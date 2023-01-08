# Computational prediction of the phenotypic effect of flavonoids on adiponectin biosynthesis
The scripts and source data used to build RAS model to predict the adiponectin-secretion-promoting phenotype.

<center><img src="./data/Graphical 400x.png" width="100%" height="100%"></center>

### Overview
This repository includes the scripts and source data for
* Collecting information on phytochemicals
* Extracting biological activity results
* Determining biological test outcomes (the decision algorithm)
* Aggregating chemical similarity and BTOs
* Training random forest-based nuclear receptor activity classifier
* Building single-cell transcriptomics-based linear regression model

### Application
1. Required packages
* RandomForest
* dplyr
* fingerprint
* rcdk
2. Load source
``` source("https://raw.githubusercontent.com/ann081993/ras/master/model/ras.R") ```
3. Functions
1) RASAR matrix
``` rasar_mat <- rasar(smiles) ```
2) RAS model
``` predicted <- ras(rasar_mat) ```

---

Please cite:
An et al., Computational prediction of the phenotypic effect of flavonoids on adiponectin biosynthesis, manuscript in preparation
