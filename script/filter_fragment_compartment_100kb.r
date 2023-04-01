#!/usr/bin/Rscript
library("optparse")
library(tools)
#conda activate delfi2
#usage: Rscript filter_fragment_compartment_100kb.r -g hg19 -a AB.rds -l ../DELFI/DELFI_Introm.R -f work/cf/7764111756397e5bd3d1375a80d7fa/02_fragmentGC/EE88211.hg19.frag.tsv_frags.rds
#in docker test: Rscript filter_fragment_compartment_100kb.r -g hg19 -a AB.rds -l iDELFI.R -f 02_fragmentGC/EE88147_frags.rds

#A/B compartments
option_list = list(
    make_option(c("-a", "--ab_file"), type="character", default="AB.rds", 
              help="A/B compartment file [default= %default]", metavar="character"),
    make_option(c("-g", "--hg_version"), type="character", default="hg19", 
              help="human genome reference version [default= %default]", metavar="character"),  
    make_option(c("-l", "--lib"), type="character", default="DELFI/DELFI_Introm.R", 
              help="path to function library", metavar="character"),                   
    make_option(c("-f", "--file"), type="character", default="hg19", 
              help="input file after fragmentGC step [default= %default]", metavar="character"), 
    make_option(c("-b", "--filter"), type="character", default="filters.hg19.rda", 
              help="filter file from AB step [default= %default]", metavar="character")                                                              
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

AB = readRDS(opt$ab)
i=opt$file
path_input_frag="./02_fragmentGC"
source(opt$lib) # # Use exon regions

#03 - 2 - make 100kb bin for each sample 
id=strsplit(basename(file_path_sans_ext(i))[1],"[.]")[[1]][1]
fragfile= i
fragments=NULL
hg_version=opt$hg_version
outdir="./03_count_100kb"
####
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
        #load("./filters.hg19.rda")
        load(opt$filter)
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
}
cat("# DELFI : filter_fragment_compartment_100kb : file ",id," Done!\n")

