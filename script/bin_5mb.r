#!/usr/bin/Rscript

#06 - combine 100kb bin into 5 mb bin

library("optparse")
library(tools)
#conda activate delfi2

#usage: Rscript bin_5mb.r -l ../DELFI/DELFI_Introm.R -f 05_combine_100kb_bin/bins_100kbcompartments.rds

option_list = list(
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),                   
    make_option(c("-f", "--file"), type="character", default="bins_100kbcompartments.rds", 
              help="input file after 100 kb combination step [default= %default]", metavar="character")                                                         
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source(opt$lib) #

dir.create("06_combine_5mb_bin")
bin_5mb(bin100kb=opt$file,outfile="06_combine_5mb_bin/bins_5mbcompartments.rds")

