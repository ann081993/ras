# Relative adiponectin score (RAS) model
The scripts and source data used to build RAS model to predict the adiponectin-secretion-promoting phenotype, described in [An et al., 2023 JCIM](https://doi.org/10.1021/acs.jcim.3c00033).

<center><img src="./data/GA.png" width="80%" height="80%"></center>

### Motivation
As a phenotypic assay for screening drug candidates against human metabolic diseases, the adipogenesis model of human bone marrow mesenchymal stem cells was used by measuring *adiponectin biosynthesis*. Notably, the relative downregulation of adiponectin, cytokine produced by adipocytes, and its cellular signalling in various metabolic diseases has been previously reported, such as in atherosclerosis, type II diabetes, and nonalcoholic fatty liver disease. Hypoadiponectinemia has also been reported in human cancers such as colorectal and breast cancers. Thus, compound-induced upregulation of adiponectin secretion has been suggested as a promising therapeutic approach for human metabolic diseases and cancers.    
Here, we constructed the relative adiponectin score (RAS) model to predict the adiponectin-secretion-promoting activities of flavonoid-associated phytochemicals.

### The RAS model
The relative adiponectin score (RAS) was obtained by sequentially applying two models: an nuclear receptor (NR) activity classifier (Model 1) and a single-cell transcriptomics-based multiple linear regression model (Model 2). The RAS was defined as a sum of predicted NR activity probabilities weighted by the relative NR contribution level determined in the Model 2.

### Application
1. Required packages
    + RandomForest
    + dplyr
    + fingerprint
    + rcdk
2. Loading source data   
```R
source("https://raw.githubusercontent.com/ann081993/ras/master/model/ras.R")
```
3. Running the RAS model   
    + RASAR matrix:  ```rasar_mat <- rasar(smiles)```
        * RASAR approach can be implemented by ```rasar``` function.
        * Structures of S. baicalensis-derived flavonoids are loaded as ```smiles``` in SMILES format.

    + RAS model: ```predicted <- ras(rasar_mat)```
        * The RASAR matrix produced above can be fed to ```ras``` function.
        * One can predict relative adiponectin score (RAS).

    + Output    
        * Output describes NR activities and RAS value of chemicals.

### Description
This repository includes the scripts and source data for
* Collecting information on phytochemicals
* Extracting biological activity results
* Determining biological test outcomes (the decision algorithm)
* Aggregating chemical similarity and BTOs
* Training random forest-based nuclear receptor activity classifier
* Building single-cell transcriptomics-based linear regression model

---

Please cite:
An et al., Computational Prediction of the Phenotypic Effect of Flavonoids on Adiponectin Biosynthesis. J. Chem. Inf. Model. 2023, 63, 3, 856–869
https://doi.org/10.1021/acs.jcim.3c00033
