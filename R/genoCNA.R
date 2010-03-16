# ----------------------------------------------------------------------
#  CNA studies
# ----------------------------------------------------------------------
genoCNA <- function(snpNames, chr, pos, LRR, BAF, pBs, sampleID, 
  Para=NULL, fixPara=FALSE, cnv.only=NULL, estimate.pi.r=TRUE, 
  estimate.pi.b=TRUE, estimate.trans.m=TRUE, outputSeg=TRUE, 
  outputSNP = 3, outputTag = sampleID, outputViterbi=FALSE, 
  Ds=c(1e10, 1e10, rep(1e8, 7)), pBs.alpha=0.001, contamination=TRUE, 
  normalGtp=NULL, geno.error=0.01, min.tp=1e-4, max.diff=0.1, 
  distThreshold=1e6, transB=c(0.5,.05,.05,0.1,0.1,.05,.05,.05,.05), 
  epsilon=0.005, K=5, maxIt=200, seg.nSNP=3, traceIt=5){

  DEBUG = FALSE
  
  if(DEBUG){
    Para=NULL; fixPara=FALSE; cnv.only=NULL; estimate.pi.r=TRUE; 
    estimate.pi.b=TRUE; estimate.trans.m=TRUE; outputSeg=TRUE; 
    outputSNP = 3; outputTag = sampleID; outputViterbi=TRUE; 
    Ds=c(1e10, 1e10, rep(1e8, 7)); pBs.alpha=0.001; contamination=TRUE; 
    normalGtp=normalGtp; geno.error=0.01; min.tp=1e-4; max.diff=0.1; 
    distThreshold=1e6; transB=c(0.5,.05,.05,0.1,0.1,.05,.05,.05,.05);
    epsilon=0.005; K=5; maxIt=200; seg.nSNP=3; traceIt=5; CNA=TRUE;
    
    files = list.files("~/research/R/genoCN/R/", "R")
    for(i in 1:length(files)){
      source(sprintf("~/research/R/genoCN/R/%s", files[i]))
    }
  }
  
  # ------------------------------------------------
  ## initialize parameters
  # ------------------------------------------------
    
  if(is.null(Para)){
    data(init.Para.CNA)
    Para = init.Para.CNA
  }

  # ------------------------------------------------
  ## sequence length
  # ------------------------------------------------
  L = length(pos)
  if(length(LRR) != L){
    stop("LRR must has the same length as pos\n")
  }
  if(length(BAF) != L){
    stop("BAF must has the same length as pos\n")
  }
  if(length(pBs) != L){
    stop("pBs must has the same length as pos\n")
  }
  if(!is.null(chr) && length(chr) != L){
    stop("chr must be null or a vector of the same length as pos\n")
  }
  
  if(is.null(cnv.only)){
    cnv.only = rep(FALSE, L)
  }else{
    if(length(cnv.only) != L){
      stop("cnv.only must be null or a vector of the same length as pos\n")
    }
    if(!is.logical(cnv.only)){
      stop("cnv.only must be a logical vector\n")
    }
  }
  
  # ------------------------------------------------
  ## check pBs, population B allele frequency
  # ------------------------------------------------
  pBs[which(is.na(pBs))] = 0.5
  pBs[which(cnv.only)] = 0.5
  
  if(any(pBs > 1 | pBs < 0)){
    stop("pBs must be between 0 and 1\n")
  }
  
  pBs[which(pBs < pBs.alpha)] = pBs.alpha
  pBs[which(pBs > 1 - pBs.alpha)] = 1 - pBs.alpha
  
  # ------------------------------------------------
  ## check chromosome and position information
  # ------------------------------------------------

  if(!is.null(chr)){
    if(!all(chr %in% c(as.character(1:90), "X"))){
      stop("invalid chromosome, valid values include 1 to 90 and X\n")
    }
    chr[chr=="X"] = "91"
    chr = as.numeric(chr)
    
    if(any(chr != sort(chr))){
      stop("data should be sorted from chromsome 1 to 90, and then chromsome X")
    }
  }

  if(any(diff(pos)[which(diff(chr)==0)] < 0)){
    stop("SNPs are not ordered by their chromatin position\n")
  }
  
  # ------------------------------------------------
  ## check normal Gtp
  # ------------------------------------------------

  if(is.null(normalGtp)){
    normalGtp=rep(-1,length(pBs))
  }
  
  if(!is.numeric(normalGtp)){
    stop("normalGtp should be a NULL or a numeric vector\n")
  }
  
  if(length(normalGtp) != length(pBs)){
    stop("length of normalGtp does not match length of p.Bs\n")
  }
  
  normalGtp[which(is.na(normalGtp))] = -1
  
  if(any(! unique(normalGtp) %in% c(-1,0,1,2))){
    stop("invalid value in normalGtp\n")
  }

  # ------------------------------------------------
  # transB
  # ------------------------------------------------
    
  Para[["transB"]] = transB

  # ------------------------------------------------
  # minimum transition probability
  # ------------------------------------------------
  N = 9
  if(length(min.tp)==1){
    min.tp = rep(min.tp, N)
  }
  if(length(min.tp) != N){
    stop(sprintf("length of min.tp is %d", length(min.tp)))
  }

  # ------------------------------------------------
  # check parameters
  # ------------------------------------------------
           
  Para = check.para(Para, CNA=TRUE, length(unique(chr)))
           
  # ------------------------------------------------
  ## estimate paramters
  # ------------------------------------------------
    
  if(fixPara){
    message("use fixed parameters...\n")
    Theta1 = list()
    touse = c("pi.r","mu.r","sd.r","pi.b","mu.b","sd.b","trans.m","trans.begin")
    touse = c(touse, "pi.rcn", "mu.rcn", "sd.rcn", "transB")
    for(tt in touse){
      Theta1[[tt]] = Para[[tt]]
    }
    Theta1[["CNA"]] = TRUE
    Theta1[["Ds"]]  = Ds
    if(estimate.pi.r){
      Theta1[["R.m"]] = max(LRR, na.rm=TRUE) - min(LRR, na.rm=TRUE)
    }else{
      Theta1[["R.m"]] = 1
    }
    Theta1[["contamination"]] = contamination
    Theta1[["geno.error"]] = geno.error
    Theta1[["pBs.alpha"]] = pBs.alpha
    class(Theta1) = "para.genoCN"
  }else{
    Theta1 = estimate.para(chr, pos, LRR, BAF, pBs, Para, cnv.only,      
        estimate.pi.r, estimate.pi.b, estimate.trans.m, pBs.alpha, 
        Ds, min.tp, max.diff, distThreshold, epsilon, K, maxIt, 
        traceIt, CNA=TRUE, contamination=contamination, 
        normalGtp=normalGtp)
  }
  
  # ------------------------------------------------
  ## outut results
  # ------------------------------------------------
    
  uchrs = sort(unique(chr))
  trans.begin = Theta1[["trans.begin"]]
  message("estimate copy number states for chromosome:\n")

  for(i in 1:length(uchrs)){
    chr1 = uchrs[i]
    wkp  = which(chr==chr1)
    
    if(chr1==91){ chr1 = "X" }else{ chr1=as.character(chr1) }
    message(chr1, " ", sep="")

    if(length(wkp) < 10){
      stop("less than 10 markers in chromosome", chr1, "\n")
    }

    if(i==1){ header=TRUE }else{ header=FALSE }
    
    eP = extract.postP(Theta1, pos[wkp], LRR[wkp], BAF[wkp], pBs[wkp], 
                       cnv.only[wkp], normalGtp=normalGtp[wkp], trans.begin[i,],
                       distThreshold)

    output(eP, snpNames[wkp], chr1, pos[wkp], sampleID, header, 
          output.loh=TRUE, CNA=TRUE, tag=outputTag, outputSeg=outputSeg, 
          outputSNP=outputSNP, seg.nSNP=seg.nSNP)

    if(outputViterbi){
      vi = viterbi(Theta1, snpNames[wkp], chr1, pos[wkp], 
                   LRR[wkp], BAF[wkp], pBs[wkp], 
                   cnv.only[wkp], sampleID=sampleID, 
                   normalGtp=normalGtp[wkp], trans.begin[i,], 
                   distThreshold)

      fname = sprintf("%s_segment_viterbi.txt", outputTag)

      if(header){
        colNames = c("chr", "start", "end", "state", "cn", "sample")
        colNames = c(colNames, c("snp1", "snp2", "n"))
        cat(colNames, file=fname, sep = "\t")
        cat("\n", file = fname, append=TRUE)
      }
      if(!is.null(vi)){
        wvi = which(vi$n >= seg.nSNP)
        
        if(length(wvi) > 0){
          vi  = vi[wvi,]
          write.table(vi, file = fname, quote = FALSE, 
                      sep = "\t", row.names = FALSE, 
                      col.names = FALSE, append=TRUE)
        }
        
      }
    }

  }
  message("\n")

  return(Theta1)
}
