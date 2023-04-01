#!/usr/bin/Rscript
library("optparse")
#conda activate delfi2
#usage: Rscript bin_chr_arm.r -l ../DELFI/DELFI_Introm.R -f 05_combine_100kb_bin/bins_100kbcompartments.rds

option_list = list(
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),                   
    make_option(c("-f", "--file"), type="character", default="bins_100kbcompartments.rds", 
              help="input file after 5 mb integration step [default= %default]", metavar="character")                                                                   
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source(opt$lib) 

dir.create("06_02_combine_chr")
bin_chr.arm(bin100kb=opt$file,outfile="06_02_combine_chr/bins_chr_compartments.rds")

