#!/usr/bin/Rscript
library("optparse")
#conda activate delfi2
#usage: Rscript build_feature.r -l ../DELFI/DELFI_Introm.R -f 06_combine_5mb_bin/bins_5mbcompartments.rds -b 06_02_combine_chr/bins_chr_compartments.rds -t 07_summarize_data/summary_data.csv.rds

option_list = list(
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),                   
    make_option(c("-f", "--file"), type="character", default="bins_5mbcompartments.rds", 
              help="input file after 5 mb integration step [default= %default]", metavar="character"),
    make_option(c("-b", "--bin_chr"), type="character", default="bins_chr_compartments.rds", 
              help="input file after bininig chr arm step [default= %default]", metavar="character"),              
    make_option(c("-t", "--table"), type="character", default="summary_data.csv.rds", 
              help="table of summary_data file [default= %default]", metavar="character")                                                                       
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source(opt$lib) 

dir.create("08_build_feature_dataframe")
build.feature.FSR_FSD.cov(bins_chr=opt$bin_chr,bins_5mb=opt$file,summary_tibble=opt$table,outfile="08_build_feature_dataframe/feature_FSR_FSD_cov.csv")
build_multiclass.feature.FSR_FSD.cov(bins_chr=opt$bin_chr,bins_5mb=opt$file,summary_tibble=opt$table,outfile="08_build_feature_dataframe/feature_multiclass_FSR_FSD_cov.csv")

