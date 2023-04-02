################################################################################
# v.0.1 : add meta_data function to summarize data, add Agilent_Human_Exon_V6_UTRs prob coverage and inverse_prob option to read_fragmentGC function
# v.0.0 : duplicate from DELFI_exome_improve.R
################################################################################

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(httr)
library(GenomicRanges)
library(rtracklayer)
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(Rsamtools)
library(Homo.sapiens)
class(Homo.sapiens)
library(devtools)
library(biovizBase)
library(RCurl)
library(tools)

library(multidplyr)
library(readxl)

library(caret)
library(pROC)

#Method to check chromosome boundary to avoid Error: 
#Error in loadFUN(x, seqname, ranges) : 
##trying to load regions beyond the boundaries of non-circular sequence "chr17"
check.fragment.boundary <- function(fragment,Hsapiens){
    for( i in paste0("chr", 1:22)){
        if( TRUE %in% (fragment[fragment[,"chr"]==i,][,"start"]==0)){
            convert_list= which((fragment[fragment[,"chr"]==i,][,"start"]==0) == TRUE)
            cat("# DELFI : check.fragment.boundary : Warning! Found some fragments out of reference boundary : at ",i,"at index",convert_list,"\n")
            cat("# DELFI : check.fragment.boundary : Probabry mapping error or use wrong reference genome (hg19 or hg38)\n")
            cat("# DELFI : check.fragment.boundary : Warning! we will convert 0 to 1 to fix this issue (not offset)\n")
            fragment[fragment[,"chr"]==i,][convert_list,"start"] = 1
        }else if (TRUE %in% (fragment[fragment[,"chr"]==i,][,"end"] > seqlengths(seqinfo(Hsapiens)[i])) ) {
            convert_list= which((fragment[fragment[,"chr"]==i,][,"end"] > seqlengths(seqinfo(Hsapiens)[i])) == TRUE)
            cat("# DELFI : check.fragment.boundary : Warning! Found some fragments out of reference boundary : at ",i,"at index",convert_list,"\n")
            cat("# DELFI : check.fragment.boundary : Probabry mapping error or use wrong reference genome (hg19 or hg38)\n")     
            cat("# DELFI : check.fragment.boundary : Try changing reference genome to hg19 or hg38\n")         
        }else{cat("# DELFI : check.fragment.boundary : ",i," : Pass\n")}
    }
    fragment
}

###00 hg19 or hg38 gaps & blacklisted regions
make_filtered_region <- function(version="hg19"){
    if (version=="hg19") {
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
        save(gaps.hg19, file="gaps.hg19.rda")
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

    }else if(version=="hg38"){
    ### hg38 gaps & blacklisted regions
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- "hg38"
    mySession <- browserSession()
    genome(mySession) <- genome
    gaps <- getTable(ucscTableQuery(mySession, track="gap"))
    gaps.hg38 <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                    gaps$chromEnd),
                    type=gaps$type)
    gaps.hg38 <- keepSeqlevels(gaps.hg38, paste0("chr", c(1:22, "X", "Y")),
                            pruning.mode="coarse")
    hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
    seqinfo(gaps.hg38) <- seqinfo(hsapiens)[seqlevels(gaps.hg38),]
    save(gaps.hg38, file="gaps.hg38.rda")
    # devtools::use_data(gaps.hg19, overwrite = TRUE)

    blacklisted.file <- httr::content(GET("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"))
    blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
                                col_names=c("seqnames", "start",
                                            "end", "name", "score"))
    blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
    filters.hg38 <- makeGRangesFromDataFrame(blacklisted.tib,
                                            keep.extra.columns=TRUE)
    filters.hg38 <- keepSeqlevels(filters.hg38, paste0("chr", c(1:22, "X", "Y")),
                            pruning.mode="coarse")
    seqinfo(filters.hg38) <- seqinfo(Hsapiens)[seqlevels(filters.hg38),]
    save(filters.hg38, file="filters.hg38.rda")
    # devtools::use_data(filters.hg38, overwrite = TRUE)
    }else {
        cat("# DELFI : make_filtered_region : Wrong hg version! Please specify human genome reference version hg19 or hg38!","\n")
    }
}

###01 Read GAlignmentPairs - load BAM file
read_GAlignmentPairs <- function(bamdir="./",id="bam",galpdir="./",mapqFilter=30){
    
    bamfile <- file.path(bamdir, id)
    indexed.bam <- gsub("$", ".bai", bamfile)

    if (!file.exists(indexed.bam)) {
        indexBam(bamfile)
    }

    param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,
                                            isSecondaryAlignment = FALSE,
                                            isUnmappedQuery = FALSE),
                        mapqFilter = mapqFilter)
    sample <- gsub(".bam", "", id)

    galpdir=sample
    dir.create(galpdir)
    galp.file <- file.path(galpdir, paste0(sample, ".rds"))
    galp <- readGAlignmentPairs(bamfile, param = param)
    saveRDS(galp, galp.file)
    galp
}

#02-fragment_gc 
## load data from source
## 1. from Object from GAlignmentPairs method
## 2. from fragment tsv file (FinaleDB)
## 3. from Rdata file containing object from GAlignmentPairs method
## 4. exome_version (hg19 only) contains 1.WGS (default) ,2. SureSelect_V7 , 3. illumina_V1, 4. Agilent_Human_Exon_V6_UTRs
read_fragmentGC <- function(source_GAlignmentPairs=NULL,source_file=NULL,outdir="./",id,hg_version="hg19",exome_version="WGS",inverse_prob="F" ){
    if(!dir.exists(outdir)){
        dir.create(outdir)
    }

    out_filename <- file.path(outdir, paste0(id, "_frags.rds"))
    if (!file.exists(out_filename)) {  
    file.create(out_filename)  
    if (!is.null(source_file)) {
        if(file_ext(source_file) == "rda"){
            galp <- readRDS(source_file)
            frags <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
                on.discordant.seqnames="drop")

        }else if(file_ext(source_file) == "tsv" || file_ext(source_file) == "gz" || file_ext(source_file) == "bgz"){
             if(file_ext(source_file) == "gz" || file_ext(source_file) == "bgz"){
                source_file=gzfile(source_file,'rt')
             }
            #Start here for checking script #AS
            fragment <- read.table(source_file,header=F,sep='\t')
            colnames(fragment) <- c("chr","start","end","score","strand")           
            chr_check=head(fragment[,"chr"],1)
            if (!"chr" %in% chr_check){
                cat("# DELFI : read_fragmentGC : Wrong chromosome column detected! (i.e. ",chr_check,") converting to chr",chr_check,"\n",sep="")
                fragment[,"chr"] = paste0("chr",fragment[,"chr"])
            }           

            cat("# DELFI : read_fragmentGC : Input file fragments =>",nrow(fragment),"\n")
            fragment <- fragment[fragment[,"score"] >= 30,]
            cat("# DELFI : read_fragmentGC : Input file fragments after quality filtering =>",nrow(fragment),"\n")
            fragment <- fragment[fragment[,"chr"] %in% paste0("chr", 1:22),]
            cat("# DELFI : read_fragmentGC : Input file fragments after quality removing chrX, chrY, and alternative contigs =>",nrow(fragment),"\n")

            if (hg_version == "hg19"){
                library(BSgenome.Hsapiens.UCSC.hg19)
                Hsapiens=BSgenome.Hsapiens.UCSC.hg19
            }else if(hg_version == "hg38"){
                library(BSgenome.Hsapiens.UCSC.hg38)
                Hsapiens=BSgenome.Hsapiens.UCSC.hg38
            }else{cat("# DELFI : bin_compartment_100kb : Wrong hg version! Please specify human genome reference version hg19 or hg38! \n")}

            fragment = check.fragment.boundary(fragment,Hsapiens)
            gc()
            frags <- GRanges(
                seqnames = Rle(fragment[,"chr"]),
                ranges = IRanges(start=fragment[,"start"], end = fragment[,"end"]),
                strand = Rle(strand(fragment[,"strand"])),
                score = fragment[,"score"])  
            gc()     
        } else{cat("# DELFI : read_fragmentGC : Wrong file format! (please use .rda or .tsv) \n")}

    } else if(!is.null(source_GAlignmentPairs)){
        frags <- granges(keepSeqlevels(source_GAlignmentPairs, paste0("chr", 1:22), pruning.mode="coarse"),
                on.discordant.seqnames="drop")        
    }else {cat("# DELFI : read_fragmentGC : Please specify one input fragment source \n")}

    ##### filter by exome capture covered regions
    #Default is SureSelect_V7
    if(exome_version=="SureSelect_V7"){
    exome_capture_v7 <- getURL('https://raw.githubusercontent.com/asangphukieo/BARMEN_Exome_capture/master/bed/SureSelect%20Clinical%20Research%20Exome%20V7/S31285117_hs_hg19/S31285117_Covered.clean.bed', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    exome <- read.table(textConnection(exome_capture_v7), header=FALSE)
    colnames(exome) = c("chr","start","end" )
    exome <- makeGRangesFromDataFrame(exome, keep.extra.columns=TRUE)
    if (inverse_prob == "T"){
        frag_exome <- frags[-queryHits(findOverlaps(frags, exome))] 
    }else{frag_exome <- frags[queryHits(findOverlaps(frags, exome))]}
    frags = frag_exome
    }else if(exome_version=="illumina_V1"){
    ##### filter by exome capture covered regions of whole_exome_agilent_1.1 version
    #gsutil cp gs://gatk-best-practices/somatic-b37/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.baits.interval_list .
    #or use https://raw.githubusercontent.com/asangphukieo/Exome_capture/main/HybSelOligos/whole_exome_illumina_coding_v1_whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list
    exome_capture_v1 <- getURL('https://raw.githubusercontent.com/asangphukieo/Exome_capture/main/whole_exome_illumina_coding_v1_whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    exome <- read.table(textConnection(exome_capture_v1), header=FALSE,comment.char="@")
    exome <- exome[,1:3]
    colnames(exome) = c("chr","start","end")
    exome <- makeGRangesFromDataFrame(exome, keep.extra.columns=TRUE)
    if (inverse_prob == "T"){
        frag_exome <- frags[-queryHits(findOverlaps(frags, exome))] 
    }else{frag_exome <- frags[queryHits(findOverlaps(frags, exome))]}
    frags = frag_exome  
    }else if(exome_version=="illumina_v1_plus_10bp_padding_minus_mito"){
    exome_capture_v1 <- getURL('https://raw.githubusercontent.com/asangphukieo/Exome_capture/main/whole_exome_illumina_coding_v1_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    exome <- read.table(textConnection(exome_capture_v1), header=FALSE,comment.char="@")
    exome <- exome[,1:3]
    colnames(exome) = c("chr","start","end")
    exome <- makeGRangesFromDataFrame(exome, keep.extra.columns=TRUE)
    if (inverse_prob == "T"){
        frag_exome <- frags[-queryHits(findOverlaps(frags, exome))] 
    }else{frag_exome <- frags[queryHits(findOverlaps(frags, exome))]}
    frags = frag_exome    
    }else if(exome_version=="Agilent_Human_Exon_V6_UTRs"){
    exome_capture_v1 <- getURL('https://raw.githubusercontent.com/asangphukieo/Exome_capture/main/Agilent_Human_Exon_V6_UTRs_Covered.bed', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    exome <- read.table(textConnection(exome_capture_v1), header=FALSE,comment.char="@")
    exome <- exome[,1:3]
    colnames(exome) = c("chr","start","end")
    exome <- makeGRangesFromDataFrame(exome, keep.extra.columns=TRUE)
    if (inverse_prob == "T"){
        frag_exome <- frags[-queryHits(findOverlaps(frags, exome))] 
    }else{frag_exome <- frags[queryHits(findOverlaps(frags, exome))]}
    frags = frag_exome    
    }
    
    ## filter outliers
    w.all <- width(frags)
    q.all <- quantile(w.all, c(0.001, 0.999))
    frags <- frags[which(w.all > q.all[1] & w.all < q.all[2])]
    gc()

    #instead of calculate GCcontent one time , which can cause memory issue with very large file (~1GB file requires ~70GB of memory)
    #we calculate GCcontent for each chromosome to reduce memory consumption
    #gcs <- GCcontent(Hsapiens, unstrand(frags))
    #frags$gc <- gcs
    for( i in paste0("chr", 1:22)){
        cat("# DELFI : read_fragmentGC : calculate GCcontent for",i,"\n")
        chr=keepSeqlevels(frags,i,pruning.mode="coarse")
        gcs_chr=GCcontent(Hsapiens, unstrand(chr))
        chr$gc <- gcs_chr
        if(i=="chr1"){
            grl <- chr
        }else{
            grl <- GRangesList(grl, chr)
            grl <- unlist(grl)
        }
        gc()
    } 
    frags <- grl    
    saveRDS(frags, out_filename )  
    gc()
    #frags  
}
}
#source_file="../FinaleDB_Validation_set/EE85761.hg19.frag.tsv.bgz"
#id="EE86240"
#out_filename = file.path("./02_fragmentGC/", paste0(id, "_frags.rds"))
#id="EE85761"

#03 - bin compartment
gc.correct <- function(coverage, bias) {
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage - coverage.pred + median(coverage)
}



#hic_compartments - Analysis of Hi-C data has shown that the genome can be divided into two compartments called A/B compartments. 
#These compartments are cell-type specific and are associated with open and closed chromatin.
#This file contains 100kb compartments of chromosome, we will bin the fragments following this compartment.
#this compartment was developed based on hg19 reference!
bin_compartment_100kb <- function(hg_version="hg19"){
    ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    AB <- read.table(textConnection(ABurl), header=TRUE)
    AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

    #select centromere and telomere parts 
    ##to remove from the analysis
    if (hg_version == "hg19"){
        library(BSgenome.Hsapiens.UCSC.hg19)
        Hsapiens=BSgenome.Hsapiens.UCSC.hg19
        load("./filters.hg19.rda")
        load("./gaps.hg19.rda")
        tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)] #####hg19
    }else if(hg_version == "hg38"){
        library(BSgenome.Hsapiens.UCSC.hg38)
        Hsapiens=BSgenome.Hsapiens.UCSC.hg38
        load("./filters.hg38.rda")
        load("./gaps.hg38.rda")
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
    AB
}


filter_fragment_compartment_100kb <- function(id=NULL,fragfile=NULL,fragments=NULL,hg_version="hg19",AB,outdir="./bins_100kb"){
    if (!is.null(fragfile)) {
        fragments <- readRDS(fragfile)
    }else if(!is.null(fragments)){
        fragments = fragments
    }else{cat("# DELFI : filter_fragment_compartment_100kb :Please specify fragments file! \n")}

    #id <- strsplit(basename(fragfile), "\\.")[[1]][1]
    if(!dir.exists(outdir)){
        dir.create(outdir)
    }
    filename <- file.path(outdir, paste0(id, "_count_100kb.rda"))

    if(!file.exists(filename)) {
        file.create(filename)
        # Filters Blacklist region
        if (hg_version == "hg19"){
            load("./filters.hg19.rda")
            fragments <- fragments[-queryHits(findOverlaps(fragments, filters.hg19))] #####hg19
        }else if(hg_version == "hg38"){
            load("./filters.hg38.rda")
            fragments <- fragments[-queryHits(findOverlaps(fragments, filters.hg38))]
        }else{cat("# DELFI : filter_fragment_compartment_100kb : Wrong hg version! Please specify human genome reference version hg19 or hg38! \n")}

        #Set cut-off of minimum and maximum fragment length (100,220)
        w.all <- width(fragments)
        fragments <- fragments[which(w.all >= 100 & w.all <= 400)] #AS: improve
        w <- width(fragments)

        #group of fragments by their length
        frag.list <- split(fragments, w)

        #count fragments in each length (100-400)
        counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
        if(min(w) > 100) {
            m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),
                        dimnames=list(rownames(counts), 100:(min(w)-1)))
            counts <- cbind(m0, counts)
        }
        save(AB,counts,w,fragments,frag.list,file=filename)
        #save(AB,counts,w,fragments,frag.list,file=filename, compress=TRUE, compression_level=9)
    }else{cat("# DELFI : filter_fragment_compartment_100kb : file output ",filename," exists! skip it!\n")}
}

detect_no_cov_region <- function(id=NULL,fragfile=NULL,report_out="no_cov.txt"){
    if (!is.null(fragfile)) {
        load(fragfile)
    }else{cat("# DELFI : detect_no_cov_region :Please specify fragfile ! \n")}

    #AP: check no coverage regions and remove it, and count fragment again
    cat("# DELFI : detect_no_cov_region : total compartment = ",nrow(counts),"\n")
    short <- rowSums(counts[,1:51]) #length = 100-150 (colnames(counts[,1:51]))
    medium <- rowSums(counts[,52:121]) #length = 151-220 (colnames(counts[,52:121]))
    long <- rowSums(counts[,122:ncol(counts)]) #length = 220- last column    
    no_cov_short=which(short==0)
    no_cov_medium=which(medium==0)
    no_cov_long=which(long==0)
    if(length(no_cov_short)>=1){
        cat("# DELFI : detect_no_cov_region : region with no coverage detected! at short fragment index ",no_cov_short,", these regions will be excluded! \n")
    #    AB <- AB[-no_cov]
    #    counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
    #    if(min(w) > 100) {
    #        m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),
    #                    dimnames=list(rownames(counts), 100:(min(w)-1)))
    #        counts <- cbind(m0, counts)
    #    }   
    #}
    }else{cat("# DELFI : detect_no_cov_region : short (100-150) fragment ",id," PASS \n")}

    if(length(no_cov_medium)>=1){
        cat("# DELFI : detect_no_cov_region : region with no coverage detected! at medium fragment index ",no_cov_medium,", these regions will be excluded! \n")
    }else{cat("# DELFI : detect_no_cov_region : medium (151-220) fragment ",id," PASS \n")}  

    if(length(no_cov_long)>=1){
        cat("# DELFI : detect_no_cov_region : region with no coverage detected! at long fragment index ",no_cov_long,", these regions will be excluded! \n")
    }else{cat("# DELFI : detect_no_cov_region : long (220-400) fragment ",id," PASS \n")}  
    cat(id,"\t",no_cov_short,"\t",no_cov_medium,"\t",no_cov_long,"\n",file=report_out,append = TRUE)  
    no_cov=c(unlist(no_cov_short),unlist(no_cov_medium),unlist(no_cov_long))
}

gc_correct_compartment_100kb <- function(id=NULL,fragfile=NULL,outdir="./bins_100kb",no_cov=NULL){
    if (!is.null(fragfile)) {
        load(fragfile)
    }else{cat("# DELFI : gc_correct_compartment_100kb :Please specify fragfile ! \n")}
    if (!is.null(no_cov)) {
        AB <- AB[-no_cov] 
        cat("# DELFI : gc_correct_compartment_100kb : compartment region ",no_cov,"were excluded!\n")
        counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
        if(min(w) > 100) {
            m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),
                        dimnames=list(rownames(counts), 100:(min(w)-1)))
            counts <- cbind(m0, counts)
        }           
    }
    #id <- strsplit(basename(fragfile), "\\.")[[1]][1]
    if(!dir.exists(outdir)){
        dir.create(outdir)
    }
    filename <- file.path(outdir, paste0(id, "_bin_100kb.rds"))
    if(!file.exists(filename)) {
        file.create(filename)
        
        #find mean of gc of each Hi-C region with 100kb bin
        olaps <- findOverlaps(fragments, AB)
        bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))
        bingc <- rep(NA, length(bin.list))
        bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))

        ### Get modes of fragment length ###
        Mode <- function(x) {
            ux <- unique(x)
            ux[which.max(tabulate(match(x, ux)))]
        }
        modes <- Mode(w)
        medians <- median(w)
        q25 <- quantile(w, 0.25)
        q75 <- quantile(w, 0.75)

        short <- rowSums(counts[,1:51])  #length = 100-150 (colnames(counts[,1:51]))
        medium <- rowSums(counts[,52:121])  #length = 151-220 (colnames(counts[,52:121]))
        if(ncol(counts) > 122){
            long <- rowSums(counts[,122:ncol(counts)])  #length = 200- 
        }else if(ncol(counts) == 122){
            long <- counts[,122]  
        }else{ 
            long <- rep(0,nrow(counts))} #in case no long fragment found (> 220bp)

        #Fragment Size Distribution (FSD) Ratio i, in range 100 - 220
        start=1  #start column
        df.bin5bp.corrected=vector("list", 24)
        for(i in 1:24){ # (100 - 220) /5 = 24 bin
            bin5bp = rowSums(counts[,start:(i*5)])
            df.bin5bp.corrected[[i]] <-data.frame(x=c(gc.correct(bin5bp, bingc))) 
            start= (i*5) +1            
        }
        #colnames(counts)[1] #range start        
        #colnames(counts)[121] #range end
        df.bin5bp.corrected = do.call(cbind, df.bin5bp.corrected)
        colnames(df.bin5bp.corrected) <-c(paste0("bin", 1:24))

        ratio <- short/medium
        ratio[which(is.na(ratio))]=0
        ratio[which(is.infinite(ratio))]=0
        short.corrected=gc.correct(short, bingc)
        medium.corrected=gc.correct(medium, bingc)        
        long.corrected=gc.correct(long, bingc)

        nfrags.corrected=gc.correct(short+medium, bingc)
        ratio.corrected=gc.correct(ratio, bingc)

        AB$short <- short
        AB$medium <- medium        
        AB$long <- long
        AB$ratio <- short/medium
        AB$nfrags <- short+medium+long
        AB$short.corrected <- short.corrected
        AB$medium.corrected <- medium.corrected        
        AB$long.corrected <- long.corrected
        AB$nfrags.corrected <- nfrags.corrected
        AB$ratio.corrected <- ratio.corrected

        for(i in 1:24){ mcols(AB)[paste0("bin", i)] = df.bin5bp.corrected[,paste0("bin", i)]}

        AB$mode <- modes
        AB$mean <- round(mean(w), 2)
        AB$median <- medians
        AB$quantile.25 <- q25
        AB$quantile.75 <- q75
        AB$frag.gc <- bingc
        for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] <- counts[,i]
        saveRDS(AB, filename)
    }else{cat("# DELFI : gc_correct_compartment_100kb : file output ",filename," exists! skip it!\n")}
}

#03.5 - combine bin compartment
combine_bin_compartment <- function(bindir="./bins_100kb",outfile="bins_100kbcompartments.rds"){
    files = list.files(path=bindir, pattern=".rds", all.files=FALSE,full.names=FALSE) 
    bins.list <- lapply(paste(bindir,files,sep="/"), readRDS)
    tib.list <- lapply(bins.list, as_tibble)
    ids=c()
    for (i in files){ ids=c(ids,strsplit(i,split="\\.")[[1]][1]) }
    names(tib.list) <- ids
    tib.list <- map2(tib.list, names(tib.list), ~ mutate(.x, id = .y)) %>%
        bind_rows() %>% dplyr::select(id, everything())

    tib.list <- tib.list %>% dplyr::select(-matches("X"))
    cat("# DELFI : combine_bin_compartment : Save tib.list in file bins_100kbcompartments.rds \n")

    tib.list$"id"=gsub("_frags_count_100kb_bin_100kb","",tib.list$"id")
    saveRDS(tib.list, outfile)
    tib.list
}

#04 - bin to 5mb
bin_5mb <- function(bin100kb="bins_100kbcompartments.rds",outfile="bins_5mbcompartments.rds"){
    df.fr <- readRDS(bin100kb)
    #master <- read_csv("sample_reference_test.csv")

    #df.fr2 <- inner_join(df.fr, master, by=c("id"="WGS ID"))
    df.fr2 = df.fr
    armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                "19p", "19q","20p","20q","21q","22q")
    df.fr2$arm <- factor(df.fr2$arm, levels=armlevels)

    ## Add column "combine"
    ## combine adjacent 100kb bins to form 5mb bins. We count starting from
    ## the telomeric end and remove the bin closest to the centromere if it is
    ## smaller than 5mb.
    df.fr2 <- df.fr2 %>% group_by(id, arm) %>%
        mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),
                            ceiling(rev((1:length(arm))/50) )))

    df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm, combine) %>%
        summarize(short2=sum(short, na.rm = TRUE),
                medium2=sum(medium, na.rm = TRUE),        
                long2=sum(long, na.rm = TRUE),
                short.corrected2=sum(short.corrected, na.rm = TRUE),
                medium.corrected2=sum(medium.corrected, na.rm = TRUE),
                long.corrected2=sum(long.corrected, na.rm = TRUE),
                hic.eigen=mean(eigen, na.rm = TRUE),
                gc=mean(C.G, na.rm = TRUE),
                ratio2=mean(ratio, na.rm = TRUE),
                ratio.corrected2=mean(ratio.corrected, na.rm = TRUE),                             
                nfrags2=sum(nfrags, na.rm = TRUE),
                nfrags.corrected2=sum(nfrags.corrected, na.rm = TRUE),               
                domain = median(as.integer(domain), na.rm = TRUE),
                short.var=var(short.corrected, na.rm = TRUE),
                long.var=var(long.corrected, na.rm = TRUE),
                nfrags.var=var(nfrags.corrected, na.rm = TRUE),
                mode_size=unique(mode, na.rm = TRUE),
                mean_size=unique(mean, na.rm = TRUE),
                median_size=unique(median, na.rm = TRUE),
                q25_size=unique(quantile.25, na.rm = TRUE),
                q75_size=unique(quantile.75, na.rm = TRUE),
                start=start[1],
                end=rev(end)[1],
                binsize = n())
    ### assign bins
    #df.fr3 <- inner_join(df.fr3, master, by=c("id"="WGS ID"))
    #df.fr3 <- df.fr3 %>% mutate(type = gsub(" Cancer|carcinoma", "", `Patient Type`, ignore.case=TRUE))

    ## remove the bin closest to the centromere if it is smaller than 5mb.
    df.fr3 <- df.fr3 %>% filter(binsize==50)
    df.fr3 <- df.fr3 %>% group_by(id) %>% mutate(bin = 1:length(id))

    cat("# DELFI : bin_5mb : Save df.fr3 in file bins_5mbcompartments.rds \n")
    saveRDS(df.fr3, outfile)

}

#04 - bin to 5mb
bin_chr.arm <- function(bin100kb="bins_100kbcompartments.rds",outfile="bins_chr_compartments.rds"){
    df.fr <- readRDS(bin100kb)
    #master <- read_csv("sample_reference_test.csv")

    #df.fr2 <- inner_join(df.fr, master, by=c("id"="WGS ID"))
    df.fr2 = df.fr
    armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                "19p", "19q","20p","20q","21q","22q")
    df.fr2$arm <- factor(df.fr2$arm, levels=armlevels)

    df.fr2 <- df.fr2 %>% group_by(id, arm) %>%
        mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/1),
                            ceiling(rev((1:length(arm))/1) )))

    df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm) %>%
        summarize(bin1=sum(bin1,na.rm=TRUE),
           bin2=sum(bin2,na.rm=TRUE),
           bin3=sum(bin3,na.rm=TRUE),
           bin4=sum(bin4,na.rm=TRUE),
           bin5=sum(bin5,na.rm=TRUE),
           bin6=sum(bin6,na.rm=TRUE),
           bin7=sum(bin7,na.rm=TRUE),
           bin8=sum(bin8,na.rm=TRUE),
           bin9=sum(bin9,na.rm=TRUE),
           bin10=sum(bin10,na.rm=TRUE),
           bin11=sum(bin11,na.rm=TRUE),
           bin12=sum(bin12,na.rm=TRUE),
           bin13=sum(bin13,na.rm=TRUE),
           bin14=sum(bin14,na.rm=TRUE),
           bin15=sum(bin15,na.rm=TRUE),
           bin16=sum(bin16,na.rm=TRUE),
           bin17=sum(bin17,na.rm=TRUE),
           bin18=sum(bin18,na.rm=TRUE),
           bin19=sum(bin19,na.rm=TRUE),
           bin20=sum(bin20,na.rm=TRUE),
           bin21=sum(bin21,na.rm=TRUE),
           bin22=sum(bin22,na.rm=TRUE),
           bin23=sum(bin23,na.rm=TRUE),
           bin24=sum(bin24,na.rm=TRUE))
    ### assign bins
    #df.fr3 <- inner_join(df.fr3, master, by=c("id"="WGS ID"))
    #df.fr3 <- df.fr3 %>% mutate(type = gsub(" Cancer|carcinoma", "", `Patient Type`, ignore.case=TRUE))

    df.fr3 <- df.fr3 %>% group_by(id) %>% mutate(bin = 1:length(id))

    cat("# DELFI : bin_chr.arm : Save df.fr3 in file bins_chr_compartments.rds \n")
    saveRDS(df.fr3, outfile)

}


#05 - Summarize data
summarize_data <- function(metadata="metadata_Cristiano_2019_update.csv",bin_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",outfile="sum.csv"){
    master <- read_csv(metadata)
    df.fr3 <- readRDS(bin_5mb)

    df.fr3 <- inner_join(df.fr3, master, by=c("id"="id.finaleDB"))

    healthy.median <- df.fr3 %>%
        group_by(bin) %>% 
        summarize(median.cov=median(nfrags2, na.rm=TRUE),
                median.short=median(short2, na.rm=TRUE),
                median.medium=median(medium2, na.rm=TRUE),
                median.long=median(long2, na.rm=TRUE),
                median.ratio=median(ratio2, na.rm=TRUE),
                median.corrected.cov=median(nfrags.corrected2, na.rm=TRUE),
                median.corrected.short=median(short.corrected2, na.rm=TRUE),
                median.corrected.medium=median(medium.corrected2, na.rm=TRUE),
                median.corrected.long=median(long.corrected2, na.rm=TRUE),
                median.corrected.ratio=median(ratio.corrected2, na.rm=TRUE),
                median.corrected.ratio2=median(short.corrected2/medium.corrected2, na.rm=TRUE))

    #change sample to id
    summary.df <- df.fr3 %>% ungroup() %>% group_by(id, `sample.disease`) %>%
        summarize(cov.cor=cor(nfrags2, healthy.median$median.cov, method="pearson", use="complete.obs"),
                short.cor=cor(short2, healthy.median$median.short, method="pearson", use="complete.obs"),
                medium.cor=cor(medium2, healthy.median$median.medium, method="pearson", use="complete.obs"),
                long.cor=cor(long2, healthy.median$median.long, method="pearson", use="complete.obs"),
                ratio.cor=cor(ratio2, healthy.median$median.ratio, method="pearson", use="complete.obs"),
                cov.corrected.cor=cor(nfrags.corrected2, healthy.median$median.corrected.cov, method="pearson", use="complete.obs"),
                short.corrected.cor=cor(short.corrected2, healthy.median$median.corrected.short, method="pearson", use="complete.obs"),
                medium.corrected.cor=cor(medium.corrected2, healthy.median$median.corrected.medium, method="pearson", use="complete.obs"),
                long.corrected.cor=cor(long.corrected2, healthy.median$median.corrected.long, method="pearson", use="complete.obs"),
                ratio.corrected.cor=cor(ratio.corrected2, healthy.median$median.corrected.ratio, method="pearson", use="complete.obs"),
                ratio2.corrected.cor=cor(short.corrected2/long.corrected2, healthy.median$median.corrected.ratio2, method="pearson", use="complete.obs"),
                nfrags = sum(nfrags2),
                mode_size=unique(mode_size),
                mean_size=unique(mean_size),
                median_size=unique(median_size),
                q25_size=unique(q25_size),
                q75_size=unique(q75_size),
                hqbases_analyzed = 100*sum(nfrags)*2,
                coverage = hqbases_analyzed/(504*5e6)
                )
                #note: 504 is the number of 5mb compartments
    summary.df <-inner_join(summary.df, master, by=c("id"="id.finaleDB","sample.disease"))
    summary.df$`type`<-  relevel(as.factor(summary.df$`sample.disease`), "Healthy")
    cat("# DELFI : bin_5mb : Save summary.df in file ",outfile,"\n")
    write.csv(summary.df,file=outfile,row.names=F,quote=F)
    saveRDS(summary.df, paste0(outfile,".rds"))
}

#05 - Summarize data
meta_data <- function(metadata="metadata_Cristiano_2019_update.csv",bin_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",outfile="sum.csv"){
    master <- read_csv(metadata)
    df.fr3 <- readRDS(bin_5mb)

    df.fr3 <- inner_join(df.fr3, master, by=c("id"="id.finaleDB"))
    #change sample to id
    summary.df <- df.fr3 %>% ungroup() %>% group_by(id, `sample.disease`) %>%
        summarize(
                nfrags = sum(nfrags2),
                mode_size=unique(mode_size),
                mean_size=unique(mean_size),
                median_size=unique(median_size),
                q25_size=unique(q25_size),
                q75_size=unique(q75_size),
                hqbases_analyzed = 100*sum(nfrags)*2,
                coverage = hqbases_analyzed/(504*5e6)
                )
                #note: 504 is the number of 5mb compartments
    summary.df <-inner_join(summary.df, master, by=c("id"="id.finaleDB","sample.disease"))
    cat("# DELFI : bin_5mb : Save summary.df in file ",outfile,"\n")
    write.csv(summary.df,file=outfile,row.names=F,quote=F)
    saveRDS(summary.df, paste0(outfile,".rds"))
}

#06 - train GBM model
build.feature.sxmxl <- function(bins_5mb="bins_5mbcompartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")

    #long feature
    #change sample to id
    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    #short feature
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    #medium feature
    features.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        spread(id, medium.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    features.sxmxl <- cbind(features.short, features.long, features.medium)
    colnames(features.sxl) <-c(paste0("short", 1:498), paste0("medium", 1:498), paste0("long", 1:498), ) 
    features.sxl$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.feature.sxmxl : Save feature table in file ",outfile,"\n")
    write.csv(features.sxmxl,file=outfile,row.names=F,quote=F)
    saveRDS(features.sxmxl, paste0(outfile,".rds"))    
    features.sxmxl
}     

#06 - train GBM model
build.feature.sxl <- function(bins_5mb="bins_5mbcompartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")

    #long feature
    #change sample to id
    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    #short feature
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    features.sxl <- cbind(features.long, features.short)
    colnames(features.sxl) <-c(paste0("long", 1:498), paste0("short", 1:498)) 
    features.sxl$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.feature.sxl : Save feature table in file ",outfile,"\n")
    write.csv(features.sxl,file=outfile,row.names=F,quote=F)
    saveRDS(features.sxl, paste0(outfile,".rds"))    
    features.sxl
}                        

#use in original script from Github
build.feature.nfragxs <- function(bins_5mb="bins_5mbcompartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")

    #Coverage feature
    #change sample to id
    features.cov <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.corrected2, id, bin) %>%
        spread(id, nfrags.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    #short feature
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    features.nfragxs <- cbind(features.cov, features.short)
    colnames(features.nfragxs) <-c(paste0("total", 1:498), paste0("short", 1:498)) 
    features.nfragxs$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    #features <- cbind(features.sl,
    #            as.matrix(summary.df %>% ungroup() %>%
    #                    dplyr::select(contains("Z Score"))))
    #features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")
    cat("# DELFI : build.feature.dataframe : Save feature table in file ",outfile,"\n")
    write.csv(features.nfragxs,file=outfile,row.names=F,quote=F)
    saveRDS(features.nfragxs, paste0(outfile,".rds"))    
    features.nfragxs
}

build.feature.ratio <- function(bins_5mb="bins_5mbcompartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    summary.df <- readRDS(summary_tibble)
    df.fr3 <- inner_join(df.fr3, summary.df, by="id")

    #Ratio feature
    #change sample to id
    features.ratio <- df.fr3  %>% ungroup() %>%
        dplyr::select(ratio.corrected2, id, bin) %>%
        spread(id, ratio.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%        
        t() %>%
        as.data.frame()

    colnames(features.ratio) <- paste0("ratio", 1:498)
    features.ratio$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.feature.ratio : Save feature table in file ",outfile,"\n")
    write.csv(features.ratio,file=outfile,row.names=F,quote=F)
    saveRDS(features.ratio, paste0(outfile,".rds"))    
    features.ratio
}


#06 - train GBM model
build.feature.sxl.nfrags_var <- function(bins_5mb="bins_5mbcompartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    bin_size=max(df.fr3$bin)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")

    #long feature
    #change sample to id
    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    #short feature
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    #short feature
    features.nfrags.var <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.var, id, bin) %>%
        spread(id, nfrags.var) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()

    feature.sxl.nfrags_var <- cbind(features.long, features.short,features.nfrags.var)
    colnames(feature.sxl.nfrags_var) <-c(paste0("long", 1:bin_size), paste0("short", 1:bin_size),paste0("var", 1:ncol(features.nfrags.var))) 
    feature.sxl.nfrags_var$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.feature.sxl.nfrags_var : Save feature table in file ",outfile,"\n")
    write.csv(feature.sxl.nfrags_var,file=outfile,row.names=F,quote=F)
    saveRDS(feature.sxl.nfrags_var, paste0(outfile,".rds"))    
    feature.sxl.nfrags_var
} 

build.feature.ratioxcov <- function(bins_5mb="bins_5mbcompartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    summary.df <- readRDS(summary_tibble)
    df.fr3 <- inner_join(df.fr3, summary.df, by="id")

    #Ratio feature
    #change sample to id
    features.ratio <- df.fr3  %>% ungroup() %>%
        dplyr::select(ratio.corrected2, id, bin) %>%
        spread(id, ratio.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%           
        t() %>%
        as.data.frame()

    cov = as.data.frame(summary.df[,c("coverage","nfrags")]) 
    rownames(cov)=summary.df$id

    features.ratioxcov = cbind(features.ratio,cov)
    colnames(features.ratioxcov) <- c(paste0("ratio", 1:498),"coverage","nfrags")
    features.ratioxcov$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.feature.ratio : Save feature table in file ",outfile,"\n")
    write.csv(features.ratioxcov,file=outfile,row.names=F,quote=F)
    saveRDS(features.ratioxcov, paste0(outfile,".rds"))    
    features.ratioxcov
}

build.feature.ratioxdeepfrag <- function(bins_5mb="bins_5mbcompartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    summary.df <- readRDS(summary_tibble)
    df.fr3 <- inner_join(df.fr3, summary.df, by="id")

    #Ratio feature
    #change sample to id
    features.ratio <- df.fr3  %>% ungroup() %>%
        dplyr::select(ratio.corrected2, id, bin) %>%
        spread(id, ratio.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%           
        t() %>%
        as.data.frame()

    cov = as.data.frame(summary.df) 
    rownames(cov)=summary.df$id

    features.ratioxcov = cbind(features.ratio,cov)
    colnames(features.ratioxcov) <- c(paste0("ratio", 1:498),colnames(summary.df))
    features.ratioxcov$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.feature.ratio : Save feature table in file ",outfile,"\n")
    write.table(features.ratioxcov,file=outfile,row.names=F,quote=F,sep="\t")
    saveRDS(features.ratioxcov, paste0(outfile,".rds"))    
    features.ratioxcov
}


#use in original script from Github
build.feature.FSR_FSD <- function(bins_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",bins_chr="06_02_combine_chr/bins_chr_compartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    bins_chr <- readRDS(bins_chr)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")
    bins_chr <- inner_join(bins_chr, summary.df, by="id")

    #Coverage feature
    #change sample to id
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.short[is.na(features.short)] <- 0

    features.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        spread(id, medium.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.medium[is.na(features.medium)] <- 0

    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.long[is.na(features.long)] <- 0

    #FSD feature 20659
    df.bin5bp = vector("list", 24)
    for(i in 1:24){ 
        #df.bin5bp[[i]] <-data.frame(x=bins_chr[paste0("bin", i)]) 
        df.bin5bp[[i]] <- bins_chr %>% ungroup() %>%
            dplyr::select(paste0("bin", i), id, bin) %>%
            spread(id, paste0("bin", i)) %>%
            dplyr::select(-bin) %>% 
            na.omit() %>%
            scale() %>%
            t() %>%
            as.data.frame()
    }
    features.FSD = do.call(cbind, df.bin5bp)
    features.FSD[is.na(features.FSD)] <- 0
    colnames(features.FSD) <-c(paste0("bin", 1:936)) #24 bin x 39 arm = 936bin

    ####
    features.FSR_FSD <- cbind(features.short,features.medium,features.long, features.FSD)
    colnames(features.FSR_FSD) <-c(paste0("short", 1:504), paste0("medium", 1:504), paste0("long", 1:504), paste0("bin", 1:936)) 
    features.FSR_FSD$type <- ifelse(summary.df$sample.disease == "Healthy", "Healthy", "Cancer")

    #features <- cbind(features.sl,
    #            as.matrix(summary.df %>% ungroup() %>%
    #                    dplyr::select(contains("Z Score"))))
    #features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")
    cat("# DELFI : build.features.FSR_FSD : Save feature table in file ",outfile,"\n")
    write.csv(features.FSR_FSD,file=outfile,row.names=F,quote=F)
    saveRDS(features.FSR_FSD, paste0(outfile,".rds"))    
    features.FSR_FSD
}


#use in original script from Github
build.feature.FSR_FSD.var <- function(bins_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",bins_chr="06_02_combine_chr/bins_chr_compartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    bins_chr <- readRDS(bins_chr)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")
    bins_chr <- inner_join(bins_chr, summary.df, by="id")

    #Coverage feature
    #change sample to id
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.short[is.na(features.short)] <- 0

    features.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        spread(id, medium.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.medium[is.na(features.medium)] <- 0

    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.long[is.na(features.long)] <- 0

    #variance coverage
    features.nfrags.var <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.var, id, bin) %>%
        spread(id, nfrags.var) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.nfrags.var[is.na(features.nfrags.var)] <- 0

    #FSD feature 20659
    df.bin5bp = vector("list", 24)
    for(i in 1:24){ 
        #df.bin5bp[[i]] <-data.frame(x=bins_chr[paste0("bin", i)]) 
        df.bin5bp[[i]] <- bins_chr %>% ungroup() %>%
            dplyr::select(paste0("bin", i), id, bin) %>%
            spread(id, paste0("bin", i)) %>%
            dplyr::select(-bin) %>% 
            #na.omit() %>%
            scale() %>%
            t() %>%
            as.data.frame()
    }
    features.FSD = do.call(cbind, df.bin5bp)
    features.FSD[is.na(features.FSD)] <- 0
    colnames(features.FSD) <-c(paste0("bin", 1:ncol(features.FSD))) #24 bin x 39 arm = 936bin

    ####
    features.FSR_FSD.var <- cbind(features.short,features.medium,features.long, features.nfrags.var, features.FSD)
    colnames(features.FSR_FSD.var) <-c(paste0("short", 1:ncol(features.short)), paste0("medium", 1:ncol(features.medium)), paste0("long", 1:ncol(features.long)),paste0("var", 1:ncol(features.nfrags.var)), paste0("bin", 1:ncol(features.FSD))) 
    features.FSR_FSD.var$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    cat("# DELFI : ",paste0("n short feature = ",ncol(features.short)),"\n")
    cat("# DELFI : ",paste0("n medium feature = ", ncol(features.medium)),"\n")
    cat("# DELFI : ",paste0("n long feature = ", ncol(features.long)),"\n")
    cat("# DELFI : ",paste0("n variance coverage feature = ", ncol(features.nfrags.var)),"\n")
    cat("# DELFI : ",paste0("n FSD feature = ", ncol(features.FSD)),"\n")

    #features <- cbind(features.sl,
    #            as.matrix(summary.df %>% ungroup() %>%
    #                    dplyr::select(contains("Z Score"))))
    #features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")
    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    write.csv(features.FSR_FSD.var,file=outfile,row.names=F,quote=F)
    saveRDS(features.FSR_FSD.var, paste0(outfile,".rds"))    
    features.FSR_FSD.var
}

#bins_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds"
#bins_chr="06_02_combine_chr/bins_chr_compartments.rds"
#summary_tibble="07_summarize_data/summary_data.csv.rds"
#use in original script from Github
build.feature.FSR_FSD.cov <- function(bins_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",bins_chr="06_02_combine_chr/bins_chr_compartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    bins_chr <- readRDS(bins_chr)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")
    bins_chr <- inner_join(bins_chr, summary.df, by="id")

    #Coverage feature
    #change sample to id
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.short[is.na(features.short)] <- 0

    features.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        spread(id, medium.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.medium[is.na(features.medium)] <- 0

    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.long[is.na(features.long)] <- 0

    features.cov <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.corrected2, id, bin) %>%
        spread(id, nfrags.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.cov[is.na(features.cov)] <- 0

    #FSD feature 20659
    df.bin5bp = vector("list", 24)
    for(i in 1:24){ 
        #df.bin5bp[[i]] <-data.frame(x=bins_chr[paste0("bin", i)]) 
        df.bin5bp[[i]] <- bins_chr %>% ungroup() %>%
            dplyr::select(paste0("bin", i), id, bin) %>%
            spread(id, paste0("bin", i)) %>%
            dplyr::select(-bin) %>% 
            #na.omit() %>%
            scale() %>%
            t() %>%
            as.data.frame()
    }
    features.FSD = do.call(cbind, df.bin5bp)
    features.FSD[is.na(features.FSD)] <- 0
    colnames(features.FSD) <-c(paste0("bin", 1:ncol(features.FSD))) #24 bin x 39 arm = 936bin

    ####
    features.FSR_FSD.var <- cbind(features.short,features.medium,features.long, features.cov, features.FSD)
    colnames(features.FSR_FSD.var) <-c(paste0("short", 1:ncol(features.short)), paste0("medium", 1:ncol(features.medium)), paste0("long", 1:ncol(features.long)),paste0("cov", 1:ncol(features.cov)), paste0("bin", 1:ncol(features.FSD))) 
    features.FSR_FSD.var$type <- ifelse(summary.df$sample.disease == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    cat("# DELFI : ",paste0("n short feature = ",ncol(features.short)),"\n")
    cat("# DELFI : ",paste0("n medium feature = ", ncol(features.medium)),"\n")
    cat("# DELFI : ",paste0("n long feature = ", ncol(features.long)),"\n")
    cat("# DELFI : ",paste0("n coverage feature = ", ncol(features.cov)),"\n")    
    cat("# DELFI : ",paste0("n FSD feature = ", ncol(features.FSD)),"\n")
    cat("# DELFI : n total features = ",ncol(features.FSR_FSD.var),"\n")
    #features <- cbind(features.sl,
    #            as.matrix(summary.df %>% ungroup() %>%
    #                    dplyr::select(contains("Z Score"))))
    #features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")
    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    write.csv(features.FSR_FSD.var,file=outfile,row.names=F,quote=F)
    saveRDS(features.FSR_FSD.var, paste0(outfile,".rds"))    
    features.FSR_FSD.var
}


#use in original script from Github
build_multiclass.feature.FSR_FSD.cov <- function(bins_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",bins_chr="06_02_combine_chr/bins_chr_compartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    bins_chr <- readRDS(bins_chr)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")
    bins_chr <- inner_join(bins_chr, summary.df, by="id")

    #Coverage feature
    #change sample to id
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.short[is.na(features.short)] <- 0

    features.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        spread(id, medium.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.medium[is.na(features.medium)] <- 0

    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.long[is.na(features.long)] <- 0

    features.cov <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.corrected2, id, bin) %>%
        spread(id, nfrags.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.cov[is.na(features.cov)] <- 0

    #FSD feature 20659
    df.bin5bp = vector("list", 24)
    for(i in 1:24){ 
        #df.bin5bp[[i]] <-data.frame(x=bins_chr[paste0("bin", i)]) 
        df.bin5bp[[i]] <- bins_chr %>% ungroup() %>%
            dplyr::select(paste0("bin", i), id, bin) %>%
            spread(id, paste0("bin", i)) %>%
            dplyr::select(-bin) %>% 
            #na.omit() %>%
            scale() %>%
            t() %>%
            as.data.frame()
    }
    features.FSD = do.call(cbind, df.bin5bp)
    features.FSD[is.na(features.FSD)] <- 0
    colnames(features.FSD) <-c(paste0("bin", 1:ncol(features.FSD))) #24 bin x 39 arm = 936bin

    ####
    features.FSR_FSD.var <- cbind(features.short,features.medium,features.long, features.cov, features.FSD)
    colnames(features.FSR_FSD.var) <-c(paste0("short", 1:ncol(features.short)), paste0("medium", 1:ncol(features.medium)), paste0("long", 1:ncol(features.long)),paste0("cov", 1:ncol(features.cov)), paste0("bin", 1:ncol(features.FSD))) 
    features.FSR_FSD.var$type <- summary.df$sample.disease # ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    cat("# DELFI : ",paste0("n short feature = ",ncol(features.short)),"\n")
    cat("# DELFI : ",paste0("n medium feature = ", ncol(features.medium)),"\n")
    cat("# DELFI : ",paste0("n long feature = ", ncol(features.long)),"\n")
    cat("# DELFI : ",paste0("n coverage feature = ", ncol(features.cov)),"\n")    
    cat("# DELFI : ",paste0("n FSD feature = ", ncol(features.FSD)),"\n")
    cat("# DELFI : n total features = ",ncol(features.FSR_FSD.var),"\n")
    #features <- cbind(features.sl,
    #            as.matrix(summary.df %>% ungroup() %>%
    #                    dplyr::select(contains("Z Score"))))
    #features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")
    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    write.csv(features.FSR_FSD.var,file=outfile,row.names=F,quote=F)
    saveRDS(features.FSR_FSD.var, paste0(outfile,".rds"))    
    features.FSR_FSD.var
}

#use in original script from Github
build.feature.FSR_FSD.cov.var <- function(bins_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",bins_chr="06_02_combine_chr/bins_chr_compartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    bins_chr <- readRDS(bins_chr)
    summary.df <- readRDS(summary_tibble)

    df.fr3 <- inner_join(df.fr3, summary.df, by="id")
    bins_chr <- inner_join(bins_chr, summary.df, by="id")

    #Coverage feature
    #change sample to id
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.short[is.na(features.short)] <- 0

    features.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        spread(id, medium.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.medium[is.na(features.medium)] <- 0

    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.long[is.na(features.long)] <- 0

    features.cov <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.corrected2, id, bin) %>%
        spread(id, nfrags.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.cov[is.na(features.cov)] <- 0

    #variance coverage
    features.nfrags.var <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.var, id, bin) %>%
        spread(id, nfrags.var) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.nfrags.var[is.na(features.nfrags.var)] <- 0

    #FSD feature 20659
    df.bin5bp = vector("list", 24)
    for(i in 1:24){ 
        #df.bin5bp[[i]] <-data.frame(x=bins_chr[paste0("bin", i)]) 
        df.bin5bp[[i]] <- bins_chr %>% ungroup() %>%
            dplyr::select(paste0("bin", i), id, bin) %>%
            spread(id, paste0("bin", i)) %>%
            dplyr::select(-bin) %>% 
            #na.omit() %>%
            scale() %>%
            t() %>%
            as.data.frame()
    }
    features.FSD = do.call(cbind, df.bin5bp)
    features.FSD[is.na(features.FSD)] <- 0
    colnames(features.FSD) <-c(paste0("bin", 1:ncol(features.FSD))) #24 bin x 39 arm = 936bin

    ####
    features.FSR_FSD.var <- cbind(features.short,features.medium,features.long, features.cov, features.nfrags.var, features.FSD)
    colnames(features.FSR_FSD.var) <-c(paste0("short", 1:ncol(features.short)), paste0("medium", 1:ncol(features.medium)), paste0("long", 1:ncol(features.long)),paste0("cov", 1:ncol(features.cov)), paste0("var", 1:ncol(features.nfrags.var)), paste0("bin", 1:ncol(features.FSD))) 
    features.FSR_FSD.var$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    cat("# DELFI : ",paste0("n short feature = ",ncol(features.short)),"\n")
    cat("# DELFI : ",paste0("n medium feature = ", ncol(features.medium)),"\n")
    cat("# DELFI : ",paste0("n long feature = ", ncol(features.long)),"\n")
    cat("# DELFI : ",paste0("n coverage feature = ", ncol(features.cov)),"\n")    
    cat("# DELFI : ",paste0("n variance coverage feature = ", ncol(features.nfrags.var)),"\n")
    cat("# DELFI : ",paste0("n FSD feature = ", ncol(features.FSD)),"\n")
    cat("# DELFI : n total features = ",ncol(features.FSR_FSD.var),"\n")
    #features <- cbind(features.sl,
    #            as.matrix(summary.df %>% ungroup() %>%
    #                    dplyr::select(contains("Z Score"))))
    #features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")
    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    write.csv(features.FSR_FSD.var,file=outfile,row.names=F,quote=F)
    saveRDS(features.FSR_FSD.var, paste0(outfile,".rds"))    
    features.FSR_FSD.var
}


#use in original script from Github
build.feature.FSR_FSD.cov.var2 <- function(bins_5mb="06_combine_5mb_bin/bins_5mbcompartments.rds",bins_chr="06_02_combine_chr/bins_chr_compartments.rds",summary_tibble="summary_tibble.rds",outfile="feature.csv"){
    df.fr3 <- readRDS(bins_5mb)
    bins_chr <- readRDS(bins_chr)
    summary.df <- readRDS(summary_tibble)
    df.fr3 <- inner_join(df.fr3, summary.df, by="id")
    bins_chr <- inner_join(bins_chr, summary.df, by="id")
    #Coverage feature
    #change sample to id
    features.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        spread(id, short.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.short[is.na(features.short)] <- 0

    features.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        spread(id, medium.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.medium[is.na(features.medium)] <- 0

    features.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        spread(id, long.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.long[is.na(features.long)] <- 0

    features.cov <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.corrected2, id, bin) %>%
        spread(id, nfrags.corrected2) %>%
        dplyr::select(-bin) %>% 
        #na.omit() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    features.cov[is.na(features.cov)] <- 0

    #cefficient of variation (not work because , mean can be minus)
    features.cof.cov <- df.fr3  %>% ungroup() %>%
        dplyr::select(nfrags.corrected2, id, bin) %>% 
        group_by(id) %>% summarize( cof.cov=sd(nfrags.corrected2,na.rm = FALSE)/mean(nfrags.corrected2,na.rm = FALSE)) %>%
        as.data.frame()
    features.cof.cov[is.na(features.cof.cov)] <- 0
    rownames(features.cof.cov) <- features.cof.cov$id
    features.cof.cov <- features.cof.cov[,"cof.cov",drop=F]

    features.cof.short <- df.fr3  %>% ungroup() %>%
        dplyr::select(short.corrected2, id, bin) %>%
        group_by(id) %>% summarize( cof.short=sd(short.corrected2,na.rm = FALSE)/mean(short.corrected2,na.rm = FALSE)) %>%
        as.data.frame()
    features.cof.short[is.na(features.cof.short)] <- 0
    rownames(features.cof.short) <- features.cof.short$id
    features.cof.short <- features.cof.short[,"cof.short",drop=F]

    features.cof.medium <- df.fr3  %>% ungroup() %>%
        dplyr::select(medium.corrected2, id, bin) %>%
        group_by(id) %>% summarize( cof.medium=sd(medium.corrected2,na.rm = FALSE)/mean(medium.corrected2,na.rm = FALSE)) %>%
        as.data.frame()
    features.cof.medium[is.na(features.cof.medium)] <- 0
    rownames(features.cof.medium) <- features.cof.medium$id
    features.cof.medium <- features.cof.medium[,"cof.medium",drop=F]

    features.cof.long <- df.fr3  %>% ungroup() %>%
        dplyr::select(long.corrected2, id, bin) %>%
        group_by(id) %>% summarize( cof.long=sd(long.corrected2,na.rm = FALSE)/mean(long.corrected2,na.rm = FALSE)) %>%
        as.data.frame()
    features.cof.long[is.na(features.cof.long)] <- 0
    rownames(features.cof.long) <- features.cof.long$id
    features.cof.long <- features.cof.long[,"cof.long",drop=F]

    #FSD feature 20659
    df.bin5bp = vector("list", 24)
    for(i in 1:24){ 
        #df.bin5bp[[i]] <-data.frame(x=bins_chr[paste0("bin", i)]) 
        df.bin5bp[[i]] <- bins_chr %>% ungroup() %>%
            dplyr::select(paste0("bin", i), id, bin) %>%
            spread(id, paste0("bin", i)) %>%
            dplyr::select(-bin) %>% 
            #na.omit() %>%
            scale() %>%
            t() %>%
            as.data.frame()
    }
    features.FSD = do.call(cbind, df.bin5bp)
    features.FSD[is.na(features.FSD)] <- 0
    colnames(features.FSD) <-c(paste0("bin", 1:ncol(features.FSD))) #24 bin x 39 arm = 936bin

    ####
    features.FSR_FSD.var <- cbind(features.short,features.medium,features.long, features.cov, features.cof.cov, features.cof.short, features.cof.medium, features.cof.long, features.FSD)
    colnames(features.FSR_FSD.var) <-c(paste0("short", 1:ncol(features.short)), paste0("medium", 1:ncol(features.medium)), paste0("long", 1:ncol(features.long)),paste0("cov", 1:ncol(features.cov)), paste0("cof.cov", 1:ncol(features.cof.cov)), paste0("cof.short", 1:ncol(features.cof.short)),paste0("cof.medium", 1:ncol(features.cof.medium)),paste0("cof.long", 1:ncol(features.cof.long)), paste0("bin", 1:ncol(features.FSD))) 
    features.FSR_FSD.var$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    cat("# DELFI : ",paste0("n short feature = ",ncol(features.short)),"\n")
    cat("# DELFI : ",paste0("n medium feature = ", ncol(features.medium)),"\n")
    cat("# DELFI : ",paste0("n long feature = ", ncol(features.long)),"\n")
    cat("# DELFI : ",paste0("n coverage feature = ", ncol(features.cov)),"\n")    
    cat("# DELFI : ",paste0("n coefficient of variation of coverage = ", ncol(features.cof.cov)),"\n")
    cat("# DELFI : ",paste0("n coefficient of variation of short fragment = ", ncol(features.cof.short)),"\n")
    cat("# DELFI : ",paste0("n coefficient of variation of medium fragment= ", ncol(features.cof.medium)),"\n")
    cat("# DELFI : ",paste0("n coefficient of variation of long fragment= ", ncol(features.cof.long)),"\n")
    cat("# DELFI : ",paste0("n FSD feature = ", ncol(features.FSD)),"\n")
    cat("# DELFI : n total features = ",ncol(features.FSR_FSD.var),"\n")
    #features <- cbind(features.sl,
    #            as.matrix(summary.df %>% ungroup() %>%
    #                    dplyr::select(contains("Z Score"))))
    #features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")
    cat("# DELFI : build.features.FSR_FSD.var : Save feature table in file ",outfile,"\n")
    write.csv(features.FSR_FSD.var,file=outfile,row.names=F,quote=F)
    saveRDS(features.FSR_FSD.var, paste0(outfile,".rds"))    
    features.FSR_FSD.var
}



#score :http://rstudio-pubs-static.s3.amazonaws.com/370944_96c386c03ac54ef3bec4535d49e92890.html
get_accuracy <- function(confusion_table){
  
  TP = confusion_table[2,2]
  TN = confusion_table[1,1]
  FN = confusion_table[1,2]
  FP = confusion_table[2,1]
  accuracy = round((TP + TN) / sum(TP,FP,TN,FN), 2)
  return(accuracy)
}

get_classification_error_rate <- function(confusion_table){
  
  TP = confusion_table[2,2]
  TN = confusion_table[1,1]
  FN = confusion_table[1,2]
  FP = confusion_table[2,1]
  classification_error_rate = round((FP + FN) / sum(TP,FP,TN,FN),2)
  return(classification_error_rate)
}

get_precision <- function(confusion_table){
  
  TP = confusion_table[2,2]
  TN = confusion_table[1,1]
  FN = confusion_table[1,2]
  FP = confusion_table[2,1]
  precision = round(TP / (TP + FP), 2)
  return(precision)
}

get_sensitivity <- function(confusion_table){
  
  TP = confusion_table[2,2]
  TN = confusion_table[1,1]
  FN = confusion_table[1,2]
  FP = confusion_table[2,1]
  sensitivity = round(TP / (TP + FN), 2)
  return(sensitivity)
}

get_specificity <- function(confusion_table){
  
  TP = confusion_table[2,2]
  TN = confusion_table[1,1]
  FN = confusion_table[1,2]
  FP = confusion_table[2,1]
  specificity = round(TN / (TN + FP), 2)
  return(specificity)
}

get_f1_score <- function(confusion_table){
  
  TP = confusion_table[2,2]
  TN = confusion_table[1,1]
  FN = confusion_table[1,2]
  FP = confusion_table[2,1]
  
  precision = round(TP / (TP + FP), 2)
  sensitivity = round(TP / (TP + FN), 2)
  f1_score = round((2 * precision * sensitivity) / (precision + sensitivity), 2)
  return(f1_score)
}


get_roc <- function(pred.tbl,fileout="auc.pdf"){
  df=pred.tbl
  # Define threshold values between 0 and 1, incrementing by 0.01
  threshold <- seq(0,1,0.01)
  
  Sensitivity<- c()
  Specificity<- c()
  
  # For every threshold value, determine 
  for (t in threshold){
    Sensitivity <- append(Sensitivity, sum((df$`Cancer` >= t & df$`obs` == "Cancer")) / sum(df$`obs` == "Cancer"))
    Specificity <- append(Specificity, sum((df$`Cancer` >= t & df$`obs` == "Healthy")) / sum(df$`obs` == "Healthy"))
  }
  # Push the resulted vectors to dataframe for plotting
  FPR=1-Specificity
  tmp_df <- data.frame(Sensitivity=Sensitivity, FPR=FPR)
  # Plot
  roc_plot <- ggplot(tmp_df, aes(x=FPR, y=Sensitivity, group=1)) + 
    geom_line() + 
    geom_point() + 
    geom_abline(intercept = 1, slope = 1)+
    xlab("False Positive Rate (1-Specificity)")+
    ylab("True Positive Rate (Sensitivity)")+
    scale_x_reverse()
  
  #Area Under the Curve (AUC)
  pos = df[df$`obs` == "Cancer", ]$Cancer
  neg = df[df$`obs` == "Healthy", ]$Cancer
  #auc = mean(replicate(100000, sample(pos, size=1) > sample(neg, size=1)))
  auc = replicate(1000,mean(replicate(1000, sample(pos, size=1) > sample(neg, size=1))))
  auc_value = mean(auc)
  ci95 = quantile( auc, probs=c(.025, .975) ) 
  
  pdf(paste0(fileout,".pdf"))
   print(roc_plot)
  dev.off()

  save(Sensitivity,Specificity,auc,auc_value,ci95,file=paste0(fileout,".Rdata"))
  output=paste0("AUC = ",auc_value,'\n',"CI = ",ci95[1]," ",ci95[2])
  write.csv(output,file=paste0(fileout,".csv"),row.names=FALSE,quote=F)


  return(list(plot=roc_plot, auc=auc_value)) 
}

auc_report <- function(ntree=150,depth=5,shrinkage=0.2,n.minobsinnode=20,model="09_build_model/models_sl.rds",feature="08_build_feature_dataframe/feature_sl.csv.rds",output="10_prediction_result/auc"){
    model_sl = readRDS(model)
    features.sl = readRDS(feature)

    pred.tbl <- model_sl$pred %>% filter(n.trees==ntree, interaction.depth==depth,shrinkage==shrinkage,n.minobsinnode==n.minobsinnode) %>%
        group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
    pred.tbl$sample <- rownames(features.sl)
    #pred.tbl <- inner_join(pred.tbl, features.sl, by=c("sample"="id"))

    rocauc <- get_roc(pred.tbl,fileout=output)   
    write.csv(pred.tbl ,paste0(output,"_pred.csv"),row.names=FALSE,quote=F)
    pred.tbl
}

auc_report_glm <- function(model="09_build_model/models_sl.rds",feature="08_build_feature_dataframe/feature_sl.csv.rds",output="10_prediction_result/auc"){
    model_sl = readRDS(model)
    features.sl = readRDS(feature)

    pred.tbl <- model_sl$pred %>% group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
    pred.tbl$sample <- rownames(features.sl)
    #pred.tbl <- inner_join(pred.tbl, features.sl, by=c("sample"="id"))

    rocauc <- get_roc(pred.tbl,fileout=output)   
    write.csv(pred.tbl ,paste0(output,"_pred.csv"),row.names=FALSE,quote=F)
    pred.tbl
}
############################################################
############ models_stack_exom.FSR_FDS.cov.var
############################################################