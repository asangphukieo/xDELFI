#!/usr/bin/Rscript
library("optparse")
library(tools)
#conda activate delfi2
#usage: Rscript read_fragmentGC.r -f ../INPUT/EE88147.hg19.frag.tsv.bgz -l ../DELFI/DELFI_Introm.R -o ./02_fragmentGC -e Agilent_Human_Exon_V6_UTRs
#in docker test: Rscript read_fragmentGC.r -f ../sample_data/EE88147.hg19.frag.tsv.bgz -l iDELFI.R -o ./02_fragmentGC -e Agilent_Human_Exon_V6_UTRs

Error in function (type, msg, asError = TRUE)  : SSL connection timeout
Calls: read_fragmentGC -> getURL -> curlPerform -> <Anonymous> -> fun
Execution halted

option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="path to input files", metavar="character"),
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),              
    make_option(c("-o", "--out"), type="character", default="02_fragmentGC", 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-g", "--hg_version"), type="character", default="hg19", 
              help="human genome reference version [default= %default]", metavar="character"),
    make_option(c("-e", "--exome_capture"), type="character", default="WGS", 
              help="exome capture probe version [default= %default]", metavar="character")                         
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#path_input="/home/sangphukieo/FinaleDB_721" #29 Dec 2022
#source("/home/sangphukieo/DELFI/DELFI_Introm.R") # # Use exon regions
#exome_capture="Agilent_Human_Exon_V6_UTRs"
#tmp_out="./02_fragmentGC"
#hg="hg19"

path_input=opt$file #29 Dec 2022
source(opt$lib) # # Use exon regions
exome_capture=opt$exome_capture
tmp_out=opt$out
hg=opt$hg_version
i=opt$file

print(path_input)
id=strsplit(basename(file_path_sans_ext(i))[1],"[.]")[[1]][1]
read_fragmentGC(source_GAlignmentPairs=NULL,
    source_file=path_input,
    outdir=tmp_out,
    id=id,
    hg_version=hg,
    exome_version=exome_capture)

