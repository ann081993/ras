# Computational Prediction of the Phenotypic Effect of Flavonoids on Adiponectin Biosynthesis
The scripts and source data used to build RAS model to predict the adiponectin-secretion-promoting phenotype.

<center><img src="./data/Graphical 400x.png" width="80%" height="80%"></center>

### Motivation
As a phenotypic assay for screening drug candidates against human metabolic diseases, the adipogenesis model of human bone marrow mesenchymal stem cells was used by measuring *adiponectin biosynthesis*. Notably, the relative downregulation of adiponectin, cytokine produced by adipocytes, and its cellular signalling in various metabolic diseases has been previously reported, such as in atherosclerosis, type II diabetes, and nonalcoholic fatty liver disease. Hypoadiponectinemia has also been reported in human cancers such as colorectal and breast cancers. Thus, compound-induced upregulation of adiponectin secretion has been suggested as a promising therapeutic approach for human metabolic diseases and cancers.   Here, we constructed the relative adiponectin score (RAS) model to predict the adiponectin-secretion-promoting activities of flavonoid-associated phytochemicals.

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
  3.1. RASAR matrix   
    + RASAR approach can be implemented by ```rasar``` function.
    + Exemplar structures of S. baicalensis-derived flavonoids are in ```smiles``` object.
```R
rasar_mat <- rasar(smiles)
```
  3.2. RAS model   
    + Using the RASAR matrix producted above can be fed to ```ras``` function.
    + One can predict relative adiponectin score (RAS).
```R
predicted <- ras(rasar_mat)
```
  3.3. Output   
    + Output describes NR activities and RAS value of chemicals.

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
An et al., Computational prediction of the phenotypic effect of flavonoids on adiponectin biosynthesis, manuscript in preparation
