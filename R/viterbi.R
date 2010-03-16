# ----------------------------------------------------------------------
#  find the best path, given parameter for ONE chromosome
#
# Since this is only for ONE chromosome, trans.begin is a vector
# ----------------------------------------------------------------------

viterbi <- function(Theta1, snpNames, chr1, pos, LRR, BAF, pBs, cnv.only,  
sampleID, normalGtp, trans.begin, distThreshold){

  # ----------------------------------------------------------------
  ## viterbi is the same as extract.postP at the beginning parts
  # ----------------------------------------------------------------
    
  touse = c("CNA", "Ds", "R.m", "contamination", "geno.error", "pBs.alpha")
  wmiss = which(! touse %in% names(Theta1))
  if(length(wmiss) > 0){
    stop("Missing components", touse[wmiss] , "in Theta1\n")
  }
  
  CNA = Theta1[["CNA"]]
  Ds  = Theta1[["Ds"]]
  R.m = Theta1[["R.m"]]
  contamination = Theta1[["contamination"]]
  geno.error = Theta1[["geno.error"]]
  pBs.alpha  = Theta1[["pBs.alpha"]]
  
  touse = c("pi.r","mu.r","sd.r","pi.b","mu.b","sd.b","trans.m","trans.begin")
  touse = c(touse, "pi.rcn", "mu.rcn", "sd.rcn", "transB")
  wmiss = which(! touse %in% names(Theta1))
  if(length(wmiss) > 0){
    stop("Missing components", touse[wmiss] , "in Theta1\n")
  }
  
  transB = Theta1[["transB"]]

  pi.r = Theta1[["pi.r"]]
  mu.r = Theta1[["mu.r"]]
  sd.r = Theta1[["sd.r"]]
  
  pi.rcn = Theta1[["pi.rcn"]]
  mu.rcn = Theta1[["mu.rcn"]]
  sd.rcn = Theta1[["sd.rcn"]]
  
  pi.b = Theta1[["pi.b"]]
  mu.b = Theta1[["mu.b"]]
  sd.b = Theta1[["sd.b"]]
  trans.m = Theta1[["trans.m"]]
  
  if(class(Theta1) != "para.genoCN"){
    warning("Theta1 is not of class 'para.genoCN'\n")
  }
  
  # ------------------------------------------------
  ## remove those SNPs with missing values
  # ------------------------------------------------
  L0  = length(pos)
  
  wNA = which(is.na(pos) | is.na(LRR) | is.na(BAF))
  
  pBs[which(is.na(pBs))] = 0.5
  if(!is.null(cnv.only)){
    pBs[which(cnv.only)] = 0.5
  }else{
    cnv.only = rep(FALSE, L0)
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
  ## check pBs, population B allele frequency
  # ------------------------------------------------
  if(any(pBs >1 | pBs < 0)){
    stop("pBs must be between 0 and 1\n")
  }
  pBs[which(pBs < pBs.alpha)] = pBs.alpha
  pBs[which(pBs > 1 - pBs.alpha)] = 1 - pBs.alpha
  
  # ------------------------------------------------
  # check parameters
  # ------------------------------------------------
  if(any(sd.r <= 0)){
    stop("sd.r must be greater than 0\n")
  }
  
  # ------------------------------------------------
  ## emission probabilities
  # ------------------------------------------------
    
  L = length(LRR)
  N = nrow(mu.b)
  M = ncol(mu.b)
  dims = c(L,N,M)
  
  eLRR = matrix(1, nrow=L, ncol=N)
  eBAF = matrix(1, nrow=L, ncol=N)
  em   = matrix(0, nrow=L, ncol=N)
  
  Z = .C("emiss", as.double(LRR), as.double(BAF),  as.double(pBs), 
         as.integer(normalGtp), as.double(geno.error), 
         as.double(pi.r), as.double(mu.r), as.double(sd.r),
         as.double(pi.rcn), as.double(mu.rcn), as.double(sd.rcn),
         as.double(R.m), as.double(pi.b), as.double(t(mu.b)),
         as.double(t(sd.b)), as.integer(dims), as.integer(!cnv.only), 
         as.integer(CNA), as.integer(contamination),
         eLRR=as.double(t(eLRR)), eBAF=as.double(t(eBAF)),
         em=as.double(t(em)), PACKAGE="genoCN")
  
  if(any(is.na(Z$em))){
    stop("some emission probabilities are missing\n")
  }
  
  if(any(Z$em == Inf)){
    stop("some emission probabilities are positive infinite\n")
  }
  
  em   = matrix(Z$em,   byrow=TRUE, nrow=L, ncol=N)
  eLRR = matrix(Z$eLRR, byrow=TRUE, nrow=L, ncol=N)
  eBAF = matrix(Z$eBAF, byrow=TRUE, nrow=L, ncol=N)
  
  rm(Z)
  
  #-----------------------------------------------------------------
  # all the above code are the sams as the code for extract.postP
  #-----------------------------------------------------------------
  
  path = numeric(L)
  logP = 0

  em.inf   = which(t(em)==-Inf)
  em.noinf = length(em.inf)
  em[em==-Inf] = 0

  dims = c(L, N, em.noinf)

  Z = .C("viterbi", as.double(pos), as.double(t(em)), 
         as.integer(em.inf),as.double(t(trans.m)), 
         as.double(trans.begin), as.double(Ds), 
         as.double(distThreshold), as.double(transB), 
         as.integer(dims), path=as.integer(path), 
         logP=as.double(logP), PACKAGE="genoCN")
    
  path = Z$path
  logP = Z$logP

  ## return result
  wcnv = which(path!=1)
  
  if(length(wcnv) == 0){ return(NULL) }
  
  pos.wkp = pos[wcnv]
  snp.wkp = snpNames[wcnv]
  sat.wkp = path[wcnv]
  
  sat.bd  = which(sat.wkp[-1] != sat.wkp[-length(sat.wkp)])
  bd      = sort(union(which(diff(wcnv) > 1), sat.bd))
  starts  = c(1, bd+1)
  ends    = c(bd, length(wcnv))
  
  nn        = length(starts)
  seg.n     = ends - starts + 1
  seg.start = pos.wkp[starts]
  seg.end   = pos.wkp[ends]
  seg.sat   = sat.wkp[starts]
  seg.snp1  = snp.wkp[starts]
  seg.snp2  = snp.wkp[ends]
  
  cn.genoCN = function(states, CNA){
    cn1 = states
    cn1[which(states==1 | states==2)] = 2
    cn1[which(states==3)] = 0
    cn1[which(states==4)] = 1
    if(CNA){
      cn1[which(states %in% 5:6)] = 3
      cn1[which(states %in% 7:9)] = 4
    }else{
      cn1[which(states==5)] = 3
      cn1[which(states==6)] = 4
    }
    cn1
  }
  
  seg.cn  = cn.genoCN(seg.sat, CNA)
    
  seg = data.frame(chr=rep(chr1, length(seg.start)), 
                   start=seg.start, end=seg.end, state=seg.sat, 
                   cn=seg.cn, sample=rep(sampleID,nn), snp1=seg.snp1, 
                   snp2=seg.snp2, n=seg.n)
  
  seg
}


