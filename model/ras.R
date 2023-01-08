source("https://raw.githubusercontent.com/ann081993/fp_util/master/fp_util_src.R")
source("https://raw.githubusercontent.com/ann081993/fp_util/master/pubchem_util.R")
source("https://raw.githubusercontent.com/ann081993/fp_util/master/utils.R")
library(randomForest)

load()

sb_fvo$smi <- c("O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(O)C=C3O",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(O)C=C3OCC=C",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(O)C=C3OCC4=CC=CC=C4",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(O)C=C3OCC4=CC=C(Cl)C=C4",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(O)C=C3OCC4=CC=C(F)C=C4",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(OCC4=CC=CC=C4)C=C3O",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(OCC4=CC=CC=C4)C=C3OC",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(OCC4=CC=CC=C4)C=C3OCC5=CC=CC=C5",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(OCC4=CC=CC=C4)C=C3OC/C=C(C)/C",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(OC/C=C(C)/C)C=C3OCC4=CC=CC=C4",
                "O=C(C=C(C1=CC=CC=C1)O2)C3=C2C=C(OC/C=C(C)/C)C=C3OC/C=C(C)/C")

rasar <- function(smiles, method = "circular") {
        rasar_result <- NULL
        fp_query <- smi2fp(smiles, type = "ECFP6")
        
        for(n in 1:nrow(fp_query)) {
                rasar_line <- NULL
                for(f in fts) {
                        load(paste0("./tr_data/", f, ".", method))
                        fp <- rbind(fp_query[n, ], fp)
                        fp <- fp[!duplicated(fp), ]
                        sim <- fp2sim(fp, part = 1)
                        rasar_line <- c(rasar_line, max(sim))
                }
                rasar_result <- rbind(rasar_result, rasar_line)
        }
        colnames(rasar_result) <- fts
        rownames(rasar_result) <- NULL
        return(data.frame(rasar_result, stringsAsFactors = FALSE))
}
