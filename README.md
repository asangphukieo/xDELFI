# Nextflow Pipeline to generate DELFI and xDELFI feature
## Installation
1. Install docker (follows https://github.com/asangphukieo/CMUTEAM_NF_pipeline)
2. Install Nextflow (e.g. by anaconda: conda install -c bioconda nextflow)
3. Clone this repository https://github.com/asangphukieo/xDELFI
4. Corect paths in test.config

## Input files
1. preprocessed cfDNA WGS data
The pipeline requires preprocessed cfDNA WGS data in a sorted BED-like, tab-separated values (TSV) file.
The data preprocessing can follow finaleDB protocol (https://finaledb.gitbook.io/finaledb-documentation).
The input of the xDELFI pipeline is compatible with BED-like TSV file from finaleDB (Chr, Start, Stop, Strand). 

2. metadatafile containing id (link to the sample name) and phenotype (Cancer and Healthy)

## Output files
Output files are reported in RDS format and comma-separated format for normalized fragmentomic feature in 5 mb bin. Header of the file shows feature type. Example profile of long fragment features are showed below
<img width="638" alt="image" src="https://github.com/asangphukieo/xDELFI/assets/47389288/482ad836-f555-4e09-999f-5fc44f01945b">

The profile feature in RDS can be used to train machine learning model by caret package in R. For example:
```
library(caret)
feature=readRDS("input_feature.RDS")
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
                     summaryFunction = twoClassSummary)

model <- caret::train(type ~ .,
                     data = feature,
                     method = 'gbm',
                     tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                         shrinkage=0.1,
                                         n.minobsinnode=10),
                     preProcess = c("corr", "nzv"),
                      trControl = ctrl) 
```

## To run the pipeline
test the pipeline with example data
```
nextflow run run_delfi.nf -with-docker asangphukieo/xdelfi:v1.0.3 -c test.config -profile docker_test
```
