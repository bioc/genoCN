# ----------------------------------------------------------------------
#  CNV studies
# ----------------------------------------------------------------------
genoCNV <- function(snpNames, chr, pos, LRR, BAF, pBs, sampleID, 
  Para=NULL, fixPara=FALSE, cnv.only=NULL, estimate.pi.r=TRUE, 
  estimate.pi.b=FALSE, estimate.trans.m=FALSE, normLRR=TRUE, 
  outputSeg=TRUE, outputSNP=3, outputTag=sampleID, outputViterbi=FALSE, 
  Ds = c(1e6, 1e6, rep(1e5, 4)), 
  pBs.alpha=0.001, loh=FALSE, output.loh=FALSE, 
  min.tp=5e-5, max.diff=0.1, distThreshold=5000, 
  transB = c(0.995, 0.005*c(.01, .09, .8, .09, .01)), 
  epsilon=0.005, K=5, maxIt=200, seg.nSNP=3, traceIt=5){

  DEBUG = FALSE
  
  if(DEBUG){
    Para=NULL; fixPara=FALSE; estimate.pi.r=TRUE; 
    estimate.pi.b=FALSE; estimate.trans.m=FALSE; normLRR=TRUE;
    outputSeg = TRUE; outputSNP = 1; outputTag = sampleID; outputViterbi=FALSE; 
    Ds = c(1e6, 1e6, rep(1e5, 4));
    pBs.alpha=0.001; loh=FALSE; output.loh=FALSE; 
    min.tp=5e-5; max.diff=0.1; distThreshold=5000;
    transB = c(0.995, 0.005*c(.01, .09, .8, .09, .01));
    epsilon=0.005; K=5; maxIt=200; seg.nSNP=3; traceIt=5;
    CNA=FALSE;
      
    files = list.files("~/research/R/genoCN/R/", "R")
    for(i in 1:length(files)){
      source(sprintf("~/research/R/genoCN/R/%s", files[i]))
    }

    message("loh =", loh, "\n")
  }
  # ------------------------------------------------
  ## initialize parameters
  # ------------------------------------------------

  if(is.null(Para)){
    data(init.Para.CNV)
    Para = init.Para.CNV
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
  pBs[which(cnv.only)]   = 0.5
  
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
  # Ds and trans.begin
  # ------------------------------------------------
  Para[["transB"]] = transB

  # ------------------------------------------------
  # if we want to skip LOH state
  # simply set the parameter of LOH state something unrealistic
  # ------------------------------------------------
  
  if(!loh){
    transM = Para[["trans.m"]]
    transM[,2] = rep(0,nrow(transM))
    for(i in 1:nrow(transM)){
      transM[i,] = transM[i,]/sum(transM[i,])
    }
    Para[["trans.m"]] = transM
    
    transBegin = Para[["trans.begin"]]
    if(is.vector(transBegin)){ transBegin=matrix(transBegin, nrow=1) }
    transBegin[,2] = rep(0.0, nrow(transBegin))
    for(i in 1:nrow(transBegin)){
      transBegin[i,] = transBegin[i,]/sum(transBegin[i,])
    }
    Para[["trans.begin"]] = transBegin
    
    transB = Para[["transB"]]
    transB[2] = 0.0
    transB = transB/sum(transB)
    Para[["transB"]] = transB
  }
  
  # ----------------------------------------------------------------------
  # minimum transition probability
  # ----------------------------------------------------------------------
  N = 6
  if(length(min.tp)==1){
    min.tp = rep(min.tp, N)
  }
  if(length(min.tp) != N){
    stop(sprintf("length of min.tp is %d", length(min.tp)))
  }
  
  if(!loh){
    min.tp[2] = 0.0
  }
    
  # ----------------------------------------------------------------------
  # check parameters
  # ----------------------------------------------------------------------
    
  Para = check.para(Para, CNA=FALSE, length(unique(chr)))

  # ----------------------------------------------------------------------
  # normalize LRR
  # ----------------------------------------------------------------------

  if(normLRR){
    medLRR  = median(LRR[LRR < 2 & LRR > -2], na.rm=TRUE)
    LRR = LRR - medLRR
  }
   
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
    Theta1[["CNA"]] = FALSE
    Theta1[["Ds"]]  = Ds
    if(estimate.pi.r){
      Theta1[["R.m"]] = max(LRR, na.rm=TRUE) - min(LRR, na.rm=TRUE)
    }else{
      Theta1[["R.m"]] = 1
    }
    Theta1[["contamination"]] = FALSE
    Theta1[["geno.error"]] = 0.01 # not used for CNV studies
    Theta1[["pBs.alpha"]] = pBs.alpha
    class(Theta1) = "para.genoCN"
  }else{
    Theta1 = estimate.para(chr, pos, LRR, BAF, pBs, Para, cnv.only, 
       estimate.pi.r, estimate.pi.b, estimate.trans.m, pBs.alpha, 
       Ds, min.tp, max.diff, distThreshold, epsilon, K, maxIt, 
       traceIt, CNA=FALSE, contamination=FALSE, normalGtp=NULL)
  }

  Theta1[["sample"]] = sampleID

  # ------------------------------------------------
  ## output results
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
                       cnv.only[wkp], normalGtp=NULL, trans.begin[i,],
                       distThreshold)
    
    output(eP, snpNames[wkp], chr1, pos[wkp], sampleID, header, 
           output.loh, CNA=FALSE, tag=outputTag, outputSeg=outputSeg, 
           outputSNP=outputSNP, seg.nSNP=seg.nSNP)
    
    if(outputViterbi){
      vi = viterbi(Theta1, snpNames[wkp], chr1, pos[wkp], 
                   LRR[wkp], BAF[wkp], pBs[wkp], 
                   cnv.only[wkp], sampleID=sampleID, 
                   normalGtp=NULL, trans.begin[i,], 
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
