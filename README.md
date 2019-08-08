---
title: "README"
author: "Melike"
date: "08/08/2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Analysing GTEx data for ageing studies

## Data

1. On Aug 8, 2019 the following data are downloaded from [GTEx website](https://www.gtexportal.org/).

```{bash}
mkdir -p data/raw/
cd data/raw
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_Analysis_v7_Annotations_SubjectPhenotypesDD.xlsx
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt
mkdir expression
cd expression
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
gunzip GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
```

