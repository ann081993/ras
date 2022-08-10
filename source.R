#### Initiate ... ----
source("https://raw.githubusercontent.com/ann081993/fp_util/master/fp_util_src.R")
source("https://raw.githubusercontent.com/ann081993/fp_util/master/pubchem_util.R")
source("https://raw.githubusercontent.com/ann081993/fp_util/master/utils.R")
library(ggplot2)
library(ggrepel)
library(gplots)
library(viridis)
library(cluster)
library(ggpubr)
library(randomForest)
library(pROC)
library(MASS)
library(caret)
library(cowplot)
library(ggsci)
library(jsonlite)
library(export)

#### KEGG FAPs ----
full.htext <- readLines("https://www.genome.jp/kegg-bin/download_htext?htext=br08003.keg&format=htext&filedir=")
kegg <- NULL
for(ind in grep("^D      ", full.htext)) {
  kegg <- rbind(kegg,
                c(full.htext[tail(grep("^A<b>", full.htext[1:ind]), 1)],
                  full.htext[tail(grep("^B  ", full.htext[1:ind]), 1)],
                  full.htext[tail(grep("^C    ", full.htext[1:ind]), 1)],
                  full.htext[ind])
  )
}
kegg[, 1] <- gsub("^A<b>|</b>", "", kegg[, 1])
kegg[, 2] <- gsub("^B  ", "", kegg[, 2])
kegg[, 3] <- gsub("^C    ", "", kegg[, 3])
kegg[, 4] <- gsub("^D      ", "", kegg[, 4])

kegg <- data.frame(cbind(kegg[, -4], matrix(unlist(strsplit(kegg[, 4], split = "  ")), ncol = 2, byrow = TRUE)),
                   stringsAsFactors = FALSE)
colnames(kegg) <- c("superclass", "class", "subclass", "entry", "name")
kegg_fvo <- kegg[kegg$superclass == "Flavonoids", ]
kegg_fvo$smi <- sapply(kegg_fvo$name, name2smi)
kegg_fvo$smi_entry <- sapply(kegg_fvo$entry, name2smi)
kegg_fvo$cid <- sapply(kegg_fvo$name, name2cid)
kegg_fvo$cid_entry <- sapply(kegg_fvo$entry, name2cid)
rownames(kegg_fvo) <- NULL

glycon <- matches("[#8]-[#6]-1-[#6]-[#8]-[#6]-[#6](-[#8])-[#6]-1-[#8]",
                  parse.smiles(kegg_fvo$smi), return.matches = FALSE)
exc_ind <- kegg_fvo$class == "Complex flavonoids" |
  kegg_fvo$class %in% c("Biflavonoids and polyflavonoids", "Diels-Alder adduct of chalcone") |
  glycon |
  grepl("S", kegg_fvo$smi)

# composition
head(kegg_fvo)
kegg_fvo$subclass[kegg_fvo$subclass %in% c("Flavan 4-ols", "Flavan 3,4-diols", "Flavan 3-ols")] <- "Flavanols"
kegg_fvo$class[exc_ind] <- "Complex flavonoids"
kegg_fvo$subclass[exc_ind] <- "Complex flavonoids"
kegg_composition <- data.frame(table(kegg_fvo$subclass), stringsAsFactors = FALSE)
colnames(kegg_composition)[1] <- "subclass"
kegg_composition <- merge(x = kegg_composition, y = kegg_fvo[, 1:3], by = "subclass")
kegg_composition <- kegg_composition[!duplicated(kegg_composition), ]
kegg_composition <- kegg_composition[order(c(0, 1000, 500)[factor(kegg_composition$class)] + kegg_composition$Freq), ]
kegg_composition$class <- factor(kegg_composition$class, levels = c("Flavonoids", "Isoflavonoids", "Complex flavonoids"))

ggpie(kegg_composition, lab.pos = "out", 
      x = "Freq", fill = "class", palette = alpha(c("cyan4", "deepskyblue2", "black"), 0.5))

#### EFAC 1: Substructures from KEGG FAPs ----
substructures <- sapply(kegg_fvo$cid, get_substructure, MaxRecords = 500)

summary(sapply(substructures, length))
table(sapply(substructures, length))
for(i in which(sapply(substructures, length) == 500)) {
  print(i)
  substructures[[i]] <- get_substructure(kegg_fvo$cid[i], MaxRecords = 2000)
}

substructures_smi <- cid2smi_bulk(unique(unlist(substructures)))

heteroatom <- grepl("S|A|c|B|r|l|F|#|a|I|P|K|e|n|Z|o|i|d|U|u|L|g|R|f|T|t|V|M|b|Y|h|W|G|s|m|D|y|E", substructures_smi$IsomericSMILES)
substructures_smi <- substructures_smi[!heteroatom, ]

#### EFAC 2: Similar from KEGG FAPs ----
similar <- sapply(kegg_fvo$cid, get_similar, MaxRecords = 20000, Threshold = 90)

similar_smi <- cid2smi_bulk(unique(unlist(similar)))

heteroatom <- grepl("S|A|c|B|r|l|F|#|a|I|P|K|e|n|Z|o|i|d|U|u|L|g|R|f|T|t|V|M|b|Y|h|W|G|s|m|D|y|E", similar_smi$IsomericSMILES)
similar_smi <- similar_smi[!heteroatom, ]

#### EFAC integration ----
extended_smi <- rbind(substructures_smi, similar_smi)
extended_smi <- extended_smi[!duplicated(extended_smi), ]

heteroatom <- grepl("S|A|c|B|r|l|F|a|I|P|K|e|n|Z|o|i|d|U|u|L|g|R|f|T|t|V|M|b|Y|h|W|G|s|m|D|y|E", extended_smi$IsomericSMILES)
extended_smi <- extended_smi[!heteroatom, ]

parts <- ceiling(nrow(extended_smi) / 2000)
glycon <- NULL
for(p in 1:parts) {
  from = (2000 * (p - 1) + 1)
  to = (2000 * (p - 1) + ifelse((parts == p) & nrow(extended_smi) %% 2000 > 0, nrow(extended_smi) %% 2000, 2000))
  g <- matches("[#8]-[#6]-1-[#6]-[#8]-[#6]-[#6](-[#8])-[#6]-1-[#8]",
               parse.smiles(extended_smi$IsomericSMILES[from:to]), return.matches = FALSE)
  glycon <- c(glycon, g)
  cat("... Detecting glycons", from, "-", to, ":", p, "of", parts, "\n")
}
extended_smi <- extended_smi[!glycon, ]

parts <- ceiling(nrow(extended_smi) / 2000)
extended_prop <- NULL
for(p in 1:parts) {
  from = (2000 * (p - 1) + 1)
  to = (2000 * (p - 1) + ifelse((parts == p) & nrow(extended_smi) %% 2000 > 0, nrow(extended_smi) %% 2000, 2000))
  desc <- eval.desc(parse.smiles(extended_smi$IsomericSMILES[from:to]),
                    c('org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor',
                      'org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor'), verbose = TRUE)
  extended_prop <- rbind(extended_prop, desc)
  cat("... Evaluating descriptors", from, "-", to, ":", p, "of", parts, "\n")
}

quantile(extended_prop$MW, probs = c(0.1, 0.9))
quantile(extended_prop$XLogP, probs = c(0.1, 0.9))

filter_mw <- extended_prop$MW > 200.01 & extended_prop$MW < 577.20
filter_logp <- extended_prop$XLogP > -0.07 & extended_prop$XLogP < 6.31

extended_prop <- extended_prop[filter_mw & filter_logp, ]

head(extended_smi); dim(extended_smi)
head(extended_prop); dim(extended_prop)
extended_prop$IsomericSMILES <- rownames(extended_prop)
extended_prop <- merge(x = extended_prop,
                       y = extended_smi,
                       by.x = "IsomericSMILES",
                       by.y = "IsomericSMILES", all.y = FALSE)

extended_prop <- extended_prop[!(extended_prop$CID %in% kegg_fvo$cid), ]

extended_assay <- get_assaysummary_bulk(extended_prop$CID)
extended_assay <- extended_assay[extended_assay$Target.GeneID %in% c(5465, 5468, 5467, 2099, 2100), ]

extended_prop <- extended_prop[extended_prop$CID %in% extended_assay$CID, ]

#### Chemical space ----
extended_prop$group <- "Extended"
extended_prop$BTO <- "with_BTO"

kegg_fvo_assay <- kegg_fvo_assay[kegg_fvo_assay$X.Target.GeneID. %in% c(5465, 5468, 5467, 2099, 2100), ]
kegg_fvo_prop <- cbind(kegg_fvo[, "smi", drop = FALSE], kegg_fvo_prop, kegg_fvo[, "cid", drop = FALSE]) # https://www.biostars.org/p/218837/

kegg_fvo_prop$group <- "KEGG"
kegg_fvo_prop$BTO <- ifelse(kegg_fvo_prop$CID %in% kegg_fvo_assay$X.CID., "with_BTO", "wo_BTO")

df_desc <- rbind(kegg_fvo_prop, extended_prop)
summary(df_desc)

df_desc$group_BTO <- paste(df_desc$group, df_desc$BTO, sep = "_")

ggscatter(df_desc[nrow(df_desc):1, ], x = "MW", y = "XLogP", size = 1,
          fill = "group", color = "group_BTO", palette = "npg", shape = 21, ellipse = TRUE) + # viridis(2, option = "C", end = 0.8, alpha = 0.5)
  scale_color_manual(values = c("black", "gray"))

#### Assay classification ----
extended_assay <- extended_assay[extended_assay$Target.GeneID %in% c(5465, 5468, 5467, 2099, 2100), ]
kegg_fvo_assay <- kegg_fvo_assay[kegg_fvo_assay$X.Target.GeneID. %in% c(5465, 5468, 5467, 2099, 2100), ]
kegg_fvo_assay$X.Bioactivity.Outcome. <- gsub('"', '', kegg_fvo_assay$X.Bioactivity.Outcome.)
kegg_fvo_assay$X.Activity.Name. <- gsub('"', '', kegg_fvo_assay$X.Activity.Name.)
kegg_fvo_assay$X.Assay.Name. <- gsub('"', '', kegg_fvo_assay$X.Assay.Name.)
kegg_fvo_assay$X.Bioassay.Type. <- gsub('"', '', kegg_fvo_assay$X.Bioassay.Type.)
colnames(kegg_fvo_assay) <- colnames(extended_assay)

assays <- rbind(kegg_fvo_assay, extended_assay)

assays <- assays[assays$Target.GeneID %in% c(5465, 5468, 5467, 2099, 2100), ]
assays$Target.GeneID <- factor(assays$Target.GeneID, levels = c(5465, 5468, 5467, 2099, 2100))
levels(assays$Target.GeneID) <- c("PPARA", "PPARG", "PPARD", "ESR1", "ESR2")

assays$group <- ifelse(assays$CID %in% kegg_fvo$cid, "KEGG", "Extended")

ggbarplot(data.frame(table(assays[, c("group", "Target.GeneID")])),
          x = "Target.GeneID", y = "Freq", legend = "none",
          xlab = "Targets", ylab = "Counts", 
          fill = "group", palette = "npg", label = "Freq")

anames <- assays$Assay.Name
anames <- unique(anames)

type_binding <- grepl("binding|affinity|displacement", tolower(anames))
type_ago <- grepl("agonist|coactivator|transactivat", tolower(anames))
type_ant <- grepl("antagonist|inhibit|identify small molecule antagonist|antagonistic activity", tolower(anames))
assay_type <- data.frame(cbind(type_binding, type_ago, type_ant))

atypes <- rep("", length(anames))
anames[type_ant & type_binding & type_ago]; atypes[type_ant & type_binding & type_ago] <- "antagonist" # antagonist
anames[!type_ant & type_binding & type_ago]; atypes[!type_ant & type_binding & type_ago] <- "agonist" # agonist
anames[type_ant & !type_binding & type_ago]; atypes[type_ant & !type_binding & type_ago] <- "antagonist" # ago & antago
anames[type_ant & type_binding & !type_ago]; atypes[type_ant & type_binding & !type_ago] <- "binding" # binding
anames[!type_ant & !type_binding & type_ago]; atypes[!type_ant & !type_binding & type_ago] <- "agonist" # agonist
anames[type_ant & !type_binding & !type_ago]; atypes[type_ant & !type_binding & !type_ago] <- "antagonist" # antagonist & binding
anames[!type_ant & type_binding & !type_ago]; atypes[!type_ant & type_binding & !type_ago] <- "binding" # binding
anames[!type_ant & !type_binding & !type_ago]; atypes[!type_ant & !type_binding & !type_ago] <- "others" # others

atypes[grepl("electivity", anames)] <- "others"
atypes[grepl("identify small molecule antagonist", anames)] <- "antagonist"
atypes[grepl("identify small molecule agonist", anames)] <- "agonist"
atypes[grepl("Transactivation of human ", anames)] <- "agonist"
atypes[grepl("radio ligand", anames)] <- "binding"
atypes[grepl("Estrogenic", anames)] <- "agonist"

others <- anames[atypes == "others"]
atypes[anames %in% others[c(1,2,13,15,16,17,18,19,30,31,39,40)]] <- "agonist"
atypes[anames %in% others[12]] <- "binding"
atypes[anames %in% others[c(34,44,46)]] <- "antagonist"

assay_annotation <- data.frame(Assay.Name = anames, 
                               Assay.Type = atypes, 
                               stringsAsFactors = FALSE)
assays <- merge(x = assays, y = assay_annotation, by = "Assay.Name")

table(!is.na(assays$PubMed.ID))
assays$Assay.Type2 <- gsub("agonist|antagonist", "modulatory", assays$Assay.Type)
assay_composition <- data.frame(table(data.frame(cbind(!is.na(assays$PubMed.ID), assays$Assay.Type2))))
colnames(assay_composition) <- c("Literature", "Type", "Freq")
assay_composition$Type <- factor(assay_composition$Type, levels = c("modulatory", "binding", "others"))

ggpie(assay_composition[assay_composition$Literature == TRUE, ], lab.pos = "in", 
      x = "Freq", fill = "Type", palette = c("white", "gray", "black")) 

ggpie(assay_composition[assay_composition$Literature == FALSE, ], lab.pos = "in", legend = 0,
      x = "Freq", fill = "Type", palette = c("white", "gray", "black"))

#### Decision algorithm ----
aggregated <- NULL
for(cid in cids) {
  for(t in c("PPARA", "PPARG", "PPARD", "ESR1", "ESR2")) {
    agg <- 0
    
    target_cid <- assays[assays$CID == cid & assays$Target.GeneID == t, c("Bioactivity.Outcome", "Assay.Type", "Literature")]
    n_lit <- sum(target_cid$Literature)
    target_cid_lit <- target_cid[target_cid$Literature, ]
    target_cid_hts <- target_cid[!target_cid$Literature, ]
    
    if(n_lit > 0) {
      target_cid <- target_cid_lit
      n_ago_act <- sum(target_cid$Bioactivity.Outcome == "Active" & target_cid$Assay.Type == "agonist")
      n_ago_inact <- sum(target_cid$Bioactivity.Outcome == "Inactive" & target_cid$Assay.Type == "agonist")
      n_antago_act <- sum(target_cid$Bioactivity.Outcome == "Active" & target_cid$Assay.Type == "antagonist")
      n_antago_inact <- sum(target_cid$Bioactivity.Outcome == "Inactive" & target_cid$Assay.Type == "antagonist")
      n_bind_act <- sum(target_cid$Bioactivity.Outcome == "Active" & target_cid$Assay.Type == "binding")
      n_bind_inact <- sum(target_cid$Bioactivity.Outcome == "Inactive" & target_cid$Assay.Type == "binding")
      
      if(n_bind_inact >= n_bind_act & n_bind_inact > 0) agg <- -1
      if(n_antago_act >= n_antago_inact & n_antago_act > 0) agg <- -1
      if(n_ago_inact > 0) agg <- -1
      if(n_ago_act >= n_ago_inact & n_ago_act > 0) agg <- 1
    }
    if(agg == 0) {
      target_cid <- target_cid_hts
      n_ago_act <- sum(target_cid$Bioactivity.Outcome == "Active" & target_cid$Assay.Type == "agonist")
      n_ago_inact <- sum(target_cid$Bioactivity.Outcome == "Inactive" & target_cid$Assay.Type == "agonist")
      n_antago_act <- sum(target_cid$Bioactivity.Outcome == "Active" & target_cid$Assay.Type == "antagonist")
      n_antago_inact <- sum(target_cid$Bioactivity.Outcome == "Inactive" & target_cid$Assay.Type == "antagonist")
      n_bind_act <- sum(target_cid$Bioactivity.Outcome == "Active" & target_cid$Assay.Type == "binding")
      n_bind_inact <- sum(target_cid$Bioactivity.Outcome == "Inactive" & target_cid$Assay.Type == "binding")
      
      if(n_bind_inact >= n_bind_act & n_bind_inact > 0) agg <- -1
      if(n_antago_act >= n_antago_inact & n_antago_act > 0) agg <- -1
      if(n_ago_inact > 0) agg <- -1
      if(n_ago_act >= n_ago_inact & n_ago_act > 0) agg <- 1
    }
    
    aggregated <- c(aggregated, agg)
    print(length(aggregated))
  }
}

#### Read-across structure-activity relationship ---- 
# function rasar() generates RASAR matrix from the SMILES structure
rasar <- function(smiles, list_fp = list.files("./tr_data", pattern = ".ecfp$")) {
  fpmatrix <- smi2fp(smiles, method = "circular")
  rasar_result <- NULL
  for(x in list_fp) {
    load(file = paste0("./tr_data/", x))
    sim_mat <- fp2jacdis(rbind(fpmatrix, fp), part = nrow(fpmatrix))
    sim_mat[sim_mat == 1] <- 0 # remove self
    rasar_result <- rbind(rasar_result, do.call(pmax, data.frame(sim_mat)))
    cat(x, "\t", nrow(fp), "\n")
  }
  return(t(rasar_result))
}

#### Relative adipoectin score ----
targets <- c("PPARA", "PPARG", "PPARD", "ESR1", "ESR2")

adipo_tpm <- read.csv("./GSE134570_TPM.txt",
                      sep = "\t", stringsAsFactors = FALSE)
dim(adipo_tpm)
genes <- adipo_tpm$gene

# function exp_prob() implements meCDF considering gene expression and the regulatory mode
exp_prob <- function(gene) {
  ind <- which(genes %in% gene)
  x_comb <- unlist(dist_diff[ind, ]) # dist_comb
  x_diff <- unlist(dist_diff[ind, ]) # dist_diff
  
  denf <- ecdf(x_comb) #plot(denf)
  probs <- denf(x_diff)
  
  return(probs)
}

probsum <- NULL
for(t in targets) {
  tf_df_t <- tf_df[tf_df$X..TF == t, ]
  tf_pos <- tf_df_t$Target[tf_df_t$Type == "Activation"]
  tf_neg <- tf_df_t$Target[tf_df_t$Type == "Repression"]
  
  probs <- cbind(sapply(tf_pos, exp_prob), 1 - sapply(tf_neg, exp_prob)) ### negative weight !!
  probs <- probs[, colMeans(probs) < 0.65 & colMeans(probs) > 0.35]
  
  tf_probsum <- colSums(t(probs)) / ncol(probs)
  probsum <- rbind(probsum, ecdf(tf_probsum)(tf_probsum))
}

exprs_target_prob <- data.frame(t(probsum))
colnames(exprs_target_prob) <- targets

exprs_target_prob$ADIPOQ <- exp_prob("ADIPOQ")
adipo_lm <- lm(ADIPOQ ~ . + 0, exprs_target_prob) # offset = rep(0, nrow(exprs_target_prob)) (this term has no effect)
plot(exprs_target_prob$ADIPOQ, adipo_lm$fitted.values); mean(adipo_lm$residuals ^ 2); adipo_lm


