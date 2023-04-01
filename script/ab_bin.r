#!/usr/bin/Rscript

library("optparse")
library(RCurl)
library(GenomicRanges)
library(biovizBase)

#conda activate delfi2
#usage: Rscript ab_bin.r -g hg19 -l /home/apiwat/hdd/DELFI/delfi_scripts/DELFI/DELFI_Introm.R

option_list = list(
    make_option(c("-g", "--hg_version"), type="character", default="hg19", 
              help="human genome reference version [default= %default]", metavar="character"),
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character")                                               
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source(opt$lib) # # Use exon regions
hg_version=opt$hg_version

#################################################################
# make_filtered_region and bin_compartment_100kb function      ##
#################################################################


library(BSgenome.Hsapiens.UCSC.hg19)
genome <- "hg19"
mySession <- browserSession()
genome(mySession) <- genome
gaps <- getTable(ucscTableQuery(mySession, track="gap"))
gaps.hg19 <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                gaps$chromEnd),
                type=gaps$type)
gaps.hg19 <- keepSeqlevels(gaps.hg19, paste0("chr", c(1:22, "X", "Y")),
                pruning.mode="coarse")
hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
seqinfo(gaps.hg19) <- seqinfo(hsapiens)[seqlevels(gaps.hg19),]
#save(gaps.hg19, file="gaps.hg19.rda")
# devtools::use_data(gaps.hg19, overwrite = TRUE)

blacklisted.file <- httr::content(GET("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
                            col_names=c("seqnames", "start",
                                        "end", "name", "score"))
blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
filters.hg19 <- makeGRangesFromDataFrame(blacklisted.tib,
                                        keep.extra.columns=TRUE)
filters.hg19 <- keepSeqlevels(filters.hg19, paste0("chr", c(1:22, "X", "Y")),
                        pruning.mode="coarse")
seqinfo(filters.hg19) <- seqinfo(Hsapiens)[seqlevels(filters.hg19),]
save(filters.hg19, file="filters.hg19.rda")
# devtools::use_data(filters.hg19, overwrite = TRUE)
######################################################
######################################################

ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
AB <- read.table(textConnection(ABurl), header=TRUE)
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

#select centromere and telomere parts 
##to remove from the analysis
if (hg_version == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    Hsapiens=BSgenome.Hsapiens.UCSC.hg19
    #load("./filters.hg19.rda")
    #load("./gaps.hg19.rda")
    tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)] #####hg19
}else if(hg_version == "hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    Hsapiens=BSgenome.Hsapiens.UCSC.hg38
    #load("./filters.hg38.rda")
    #load("./gaps.hg38.rda")
    tcmeres <- gaps.hg38[grepl("centromere|telomere", gaps.hg38$type)]
}else{cat("# DELFI : bin_compartment_100kb : Wrong hg version! Please specify human genome reference version hg19 or hg38! \n")}

chromosomes <- GRanges(paste0("chr", 1:22),
                    IRanges(0, seqlengths(Hsapiens)[1:22]))

arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
arms <- arms[-c(25,27,29,41,43)]

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
            "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
            "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
            "19p", "19q","20p","20q","21q","22q")
arms$arm <- armlevels

#Select hic_compartments not in gap region
if (hg_version == "hg19"){
    AB <- AB[-queryHits(findOverlaps(AB, gaps.hg19))] #####hg19
}else if(hg_version == "hg38"){
    AB <- AB[-queryHits(findOverlaps(AB, gaps.hg38))] #####hg38
}    

#Select hic_compartments not in centromere|telomere region
AB <- AB[queryHits(findOverlaps(AB, arms))]
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]
seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)
AB$gc <- GCcontent(Hsapiens, AB)

## These bins had no coverage
#AB <- AB[-c(8780, 13665)]     
saveRDS(AB,"AB.rds")
cat("# DELFI : ab_bin done!\n" )