
# ----------------------------------------------------------------------
# estimate the parameters of the HMM
# ----------------------------------------------------------------------

estimate.para <- function(chr, pos, LRR, BAF, pBs, Para, cnv.only, 
  estimate.pi.r, estimate.pi.b, estimate.trans.m, pBs.alpha, Ds, 
  min.tp, max.diff, distThreshold, epsilon, K, maxIt, traceIt, CNA, 
  contamination, normalGtp=NULL, R.m=NULL, geno.error=0.01)
{
  DEBUG = FALSE
  
  if(DEBUG){
    contamination = TRUE
    geno.error = 0.01
    R.m        = NULL
    # normalGtp  = NULL
  }
  
  # ------------------------------------------------
  ## remove those SNPs with missing values
  # ------------------------------------------------
  wNA = which(is.na(pos) | is.na(LRR) | is.na(BAF))
  if(!is.null(chr)){
    wNA = union(wNA, which(is.na(chr)))
  }
    
  if(length(wNA) > 0){
    pos = pos[-wNA]
    LRR = LRR[-wNA]
    BAF = BAF[-wNA]
    pBs = pBs[-wNA]
    cnv.only = cnv.only[-wNA]

    if(!is.null(normalGtp)){
      normalGtp = normalGtp[-wNA]
    }
    
    if(!is.null(chr)){
      chr = chr[-wNA]
    }
  }
  
  if(estimate.pi.r){
    if(is.null(R.m)){
      R.m = max(LRR) - min(LRR)
    }
  }else{
    R.m = 1.0 
  }
  
  # ------------------------------------------------
  ## check normalGtp, genotype in normal tissue
  # ------------------------------------------------
  if(!CNA){
    contamination = FALSE
    normalGtp = rep(-1,length(pBs))
  }else{
    if(is.null(normalGtp)){
      normalGtp = rep(-1,length(pBs))
    }
  
    if(!is.numeric(normalGtp)){
      stop("normalGtp should be a NULL or a numeric vector\n")
    }
  
    if(length(normalGtp) != length(pBs)){
      stop("length of normalGtp does not match length of pBs\n")
    }
  
    normalGtp[which(is.na(normalGtp))] = -1
  
    if(any(! unique(normalGtp) %in% c(-1,0,1,2))){
      stop("invalid value in normalGtp\n")
    }
  }
  
  # ------------------------------------------------
  ## number of SNPs in each chromosome
  # ------------------------------------------------
  if(is.null(chr)){
    len = length(pos)
  }else{
    len = tapply(pos, chr, length) 
  }
  
  # ------------------------------------------------
  ## parameters to update
  # ------------------------------------------------

  Pats = c("pi.r", "mu.r", "sd.r", "pi.b", "mu.b", "sd.b")
  Pats = c(Pats, "pi.rcn", "mu.rcn", "sd.rcn")
  Pats = c(Pats, "trans.m", "trans.begin")

  nIt  = indicator = 0
  logL = numeric(maxIt)
  maxD = numeric(length(Pats))
  
  if(!estimate.pi.r){
    cat("fix pi.r: pi.r=(", paste(Para[["pi.r"]], collapse=", "), ")\n", sep="")
    cat("fix pi.rcn: pi.rcn=(", paste(Para[["pi.rcn"]], collapse=", "), ")\n", sep="")
  }

  if(!estimate.pi.b){
    cat("fix pi.b: pi.b=(", paste(Para[["pi.b"]], collapse=", "), ")\n", sep="")
  }

  if(!estimate.trans.m){
    cat("fix trans.m: trans.m=", "\n")
    print(Para[["trans.m"]])
    cat("\n")
  }
  
  while(indicator < K){
    nIt  = nIt + 1

    if(traceIt>0){
      if(nIt %% traceIt==0)
        cat(nIt, date(), "\n")
    }

    Theta1 = baum.welch(pos, LRR, BAF, pBs, Para, R.m, Ds, cnv.only, 
                        estimate.pi.r, estimate.pi.b, estimate.trans.m, 
                        CNA, min.tp, max.diff, len, distThreshold, 
                        contamination, normalGtp, geno.error)
        
    logL[nIt] = Theta1[["logL"]]
    for(j in 1:length(Pats)){
      Pat = Pats[j]
      maxD[j] = max(abs(Theta1[[Pat]] - Para[[Pat]]))
    }

    if(max(maxD)>epsilon){
      indicator = 0
    }else {
      indicator = indicator + 1
    }
    
    for(j in 1:length(Pats)){
      Pat = Pats[j]
      Para[[Pat]] = Theta1[[Pat]]
    }

    if(nIt >= maxIt){ break }
  }
  
  if(traceIt>0 & nIt < maxIt){
    cat("converges after", nIt, "iterations\n")
  }

  if(nIt >= maxIt){
    warning(sprintf("fail to converge after %d iterations\n", maxIt))
  }
  
  Theta1[["CNA"]] = CNA
  Theta1[["R.m"]] = R.m
  Theta1[["Ds"]]  = Ds
  Theta1[["min.tp"]] = min.tp
  Theta1[["max.diff"]] = max.diff
  Theta1[["contamination"]] = contamination
  Theta1[["geno.error"]] = geno.error
  Theta1[["epsilon"]]   = epsilon
  Theta1[["pBs.alpha"]] = pBs.alpha

  Theta1[["K"]] = K
  Theta1[["logL"]] = logL[1:nIt]

  class(Theta1) = "para.genoCN"
  Theta1
}

