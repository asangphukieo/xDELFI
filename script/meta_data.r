#!/usr/bin/Rscript

#07 - summarize with metadata

library("optparse")
#conda activate delfi2
#usage: Rscript meta_data.r -l ../DELFI/DELFI_Introm.R -f 06_combine_5mb_bin/bins_5mbcompartments.rds -t ../metadata_Cristiano_2019_update.csv
#in docker test: Rscript meta_data.r -l iDELFI.R -f 06_combine_5mb_bin/bins_5mbcompartments.rds -t ../sample_data/metadata_all_2235.csv

option_list = list(
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),                   
    make_option(c("-f", "--file"), type="character", default="bins_5mbcompartments.rds", 
              help="input file after 5 mb integration step [default= %default]", metavar="character"),
    make_option(c("-t", "--table"), type="character", default="metadata_delfi721.txt", 
              help="table of metadata file [default= %default]", metavar="character")                                                                       
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source(opt$lib) #

dir.create("07_summarize_data")
meta_data(metadata=opt$table,bin_5mb=opt$file,outfile="07_summarize_data/summary_data.csv")

