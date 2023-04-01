#!/usr/bin/Rscript

#05 - combine 100kb bin (requires memory ~220 GB for exome data of 1473 files)

library("optparse")
library(tools)
#conda activate delfi2
# 4 - script for gc correction of each 100 kb compartment
#usage: Rscript combine_bin_compartment.r -l ../DELFI/DELFI_Introm.R -d 04_bins_100kb
#in docker test: Rscript combine_bin_compartment.r -l iDELFI.R -d 04_bins_100kb

option_list = list(
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),                   
    make_option(c("-d", "--dir"), type="character", default="./04_bins_100kb", 
              help="directory of input files [default= %default]", metavar="character")                                                           
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
i=opt$file
source(opt$lib) #

dir.create("05_combine_100kb_bin")
combine_bin_compartment(bindir=opt$dir,outfile="05_combine_100kb_bin/bins_100kbcompartments.rds")
