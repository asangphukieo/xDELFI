# Nextflow Pipeline to generate DELFI and xDELFI feature

xDELFI is new machine-learning model for cancer detection using whole exome-based cfDNA profiles.<br>
This repository provides feature extraction process. Firstly, we applied three fragment size thresholds, which are short (100-150), medium (150-220), and long (>220), to count the number of fragments in each 100k bin along the genome. Secondly, overall fragment coverage is counted by summation of all fragment size in each 100k bin. Thirdly, loess regression-based approach is applied to account for GC content bias for each 100kb bin. Then, the 100k bins were combined into 5Mb bins (504 bins) for each chromosome arm. Next, the fragment size distribution (FSD) of 5bp bin in range of 100 to 220bp (24 bins) in each chromosome arm (39 arms) was calculated, and loess regression-based approach is applied to account for GC content bias for each bin. The total number of features including fragment length coverage (short, medium, long: 504x3), overall coverage (cov: 504), and FSD (bin: 24x39) is 2,952.


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
