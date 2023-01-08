library(randomForest)
library(dplyr)
source("https://raw.githubusercontent.com/ann081993/utils/master/fp_util_src.R")
source("https://raw.githubusercontent.com/ann081993/utils/master/pubchem_util.R")
source("https://raw.githubusercontent.com/ann081993/utils/master/utils.R")
download.file("https://raw.github.com/ann081993/ras/master/data/model_list.RData", destfile = "tmp")
load("tmp")
download.file("https://raw.github.com/ann081993/ras/master/data/fp_list.RData", destfile = "tmp")
load("tmp")
file.remove("tmp")


targets <- c("PPARA", "PPARG", "PPARD", "ESR1", "ESR2")
smiles <- c("C1=CC=C(C=C1)C2=CC(=O)C3=C(O2)C(=C(C=C3O)O)O",
            "COC1=C(C=C(C2=C1OC(=CC2=O)C3=CC=CC=C3)O)O",
            "COC1=C(C2=C(C(=C1)O)C(=O)C=C(O2)C3=CC=CC=C3)OC",
            "COC1=C(C2=C(C=C1O)OC(=CC2=O)C3=CC=CC=C3)O",
            "COC1=C(C(=C2C(=C1O)C(=O)C=C(O2)C3=CC=CC=C3)O)O",
            "COC1=C(C(=C2C(=C1)OC(=CC2=O)C3=CC=CC=C3)O)OC",
            "COC1=C(C(=C2C(=C1O)C(=O)C=C(O2)C3=CC=CC=C3)OC)OC",
            "COC1=C(C2=C(C(=C1)O)C(=O)C=C(O2)C3=C(C=CC=C3O)O)OC",
            "COC1=CC=CC(=C1C2=CC(=O)C3=C(C(=C(C(=C3O2)OC)OC)OC)O)O",
            "C1C(OC2=C(C1=O)C(=C(C(=C2)O)O)O)C3=CC=CC=C3",
            "COC1=C(C2=C(C=C1O)O[C@@H](CC2=O)C3=CC=CC=C3)O"
)

rasar <- function(smiles) {
        rasar_result <- NULL
        fp_query <- smi2fp(smiles, type = "ECFP6", verbose = F)
        
        for(n in 1:nrow(fp_query)) {
                rasar_line <- NULL
                for(k in 1:length(fp_list)) {
                        fp <- rbind(fp_query[n, ], fp_list[[k]])
                        fp <- fp[!duplicated(fp), ]
                        sim <- fp2sim(fp, part = 1, verbose = F)
                        rasar_line <- c(rasar_line, max(sim))
                }
                rasar_result <- rbind(rasar_result, rasar_line)
        }
        colnames(rasar_result) <- names(fp_list)
        rownames(rasar_result) <- NULL
        rasar_result <- data.frame(rasar_result, stringsAsFactors = FALSE)
        rasar_result
}

ras <- function(rasar_mat) {
        predicted <- NULL
        for(m in model_list) {
                fits <- m
                pr <- NULL
                for(n in 1:length(fits)) {
                        pr <- rbind(pr, predict(fits[[n]], rasar_mat, type = "prob")[, 2] )
                }
                predicted <- cbind(predicted, colMeans(pr))
                predicted <- data.frame(predicted)
        }
        colnames(predicted) <- paste0(names(model_list), "_pred")
        rownames(predicted) <- NULL
        predicted <- mutate(predicted, RAS = (0.21251 * PPARA_pred + 0.56570 * PPARG_pred + 0.26220 * PPARD_pred + 0.03389 * ESR1_pred - 0.08082 * ESR2_pred + 0.13819))
        predicted
}
