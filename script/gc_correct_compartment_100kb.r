#!/usr/bin/Rscript
library("optparse")
library(tools)
#conda activate delfi2
# 4 - script for gc correction of each 100 kb compartment
#usage: Rscript gc_correct_compartment_100kb.r -l ../DELFI/DELFI_Introm.R -f work/8c/c259e6818b52932af16d90031f937a/03_count_100kb/EE88211_count_100kb.rda
#in docker test: Rscript gc_correct_compartment_100kb.r -l iDELFI.R -f 03_count_100kb/EE88147_frags_count_100kb.rda

option_list = list(
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),                   
    make_option(c("-f", "--file"), type="character", default="hg19", 
              help="input file after fragmentGC step [default= %default]", metavar="character")                                                           
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
i=opt$file
source(opt$lib) #

no_cov=c(8780, 13665) #remove only two common regions that have low coverage fragments
id=strsplit(basename(file_path_sans_ext(i))[1],"[.]")[[1]][1]
filename <- file.path("./04_bins_100kb", paste0(id, "_bin_100kb.rds"))
gc_correct_compartment_100kb( id=id,fragfile=i,outdir="./04_bins_100kb",no_cov=no_cov)
cat("# DELFI : gc_correct_compartment_100kb : file output ",filename," Done!\n")

