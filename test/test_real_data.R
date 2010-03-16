library(genoCN)
setwd("~/research/CNV/data/HapMap/CEU/")

info = read.table("snp_info_hg18.txt", header=TRUE, sep="\t", as.is=TRUE)      
dat  = read.table("NA06993.txt", sep=",", header=TRUE, as.is=TRUE)

dat[1:2,]
info[1:2,]

if(any(! dat$Name %in% info$Name)){
  stop("missing SNP information \n")
}

mt  = match(info$Name, dat$Name)
wNA = which(is.na(mt))
if(length(wNA)>0){
  info = info[-wNA,]
  mt   = mt[-wNA] 
}
dat = dat[mt,]

wkp  = which(info$Chr %in% as.character(1:22))
dat  = dat[wkp,]
info = info[wkp,]

sampleID = dat$Sample[1]

Theta = genoCNV(info$Name, info$Chr, info$Position, dat$LRR, dat$BAF, 
  info$PFB, sampleID, cnv.only=(info$PFB>1), outputSeg = TRUE, 
  outputSNP = 1, outputTag = sampleID)

