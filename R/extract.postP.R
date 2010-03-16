
# ----------------------------------------------------------------------
# extract posterior probability of HMM states, copy number, or 
# genotype classes, for ONE chromosome
#
# Since this is only for ONE chromosome, trans.begin is a vector
# ----------------------------------------------------------------------

extract.postP <- function(Theta1, pos, LRR, BAF, pBs, cnv.only,  
 normalGtp, trans.begin, distThreshold){

  touse = c("CNA", "Ds", "R.m", "contamination", "geno.error", "pBs.alpha")
  wmiss = which(! touse %in% names(Theta1))
   
  if(length(wmiss) > 0){
    stop("Missing components ", touse[wmiss] , " in Theta1\n")
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
    stop("Missing components ", touse[wmiss] , " in Theta1\n")
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
    warning("Theta1 is not output of function em.update\n")
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
   
  # ------------------------------------------------
  ## forward-backward, posteriors for HMM states
  # ------------------------------------------------

  f = b = pPN = matrix(0, nrow=L, ncol=N)

  ## forward probability matrix
  f = forward(pos, em, trans.m, trans.begin, Ds, distThreshold, transB)

  ## bakward probability matrix
  b = backward(pos, em, trans.m, Ds, distThreshold, transB)

  ## posterior probability matrix
  pPN = postP(f, b)

  # ------------------------------------------------
  ## ws: the weight for each normal mixture
  # ------------------------------------------------
  
  ws = get.weights(pBs, CNA, contamination, normalGtp, geno.error)

  # ------------------------------------------------
  ## posteriors for each genotype components within 
  ## each HMM state
  # ------------------------------------------------

  if(CNA){
    H = c(3, 4, 1, 4, 4, 4, 3, 4, 4)
    idxStart = c(1,4,8,9,13,22,26,29,33)
  }else{
    H = c(3, 4, 1, 4, 4, 5)
    idxStart = c(1,4,8,9,13,17)
  }
    
  Omiga = (1:L)[which(BAF>0 & BAF<1 & (!cnv.only))]

  x = 1
  if(CNA){ MM = 36 } else{ MM = 21 }
   
  postP36N = matrix(0, L, MM)

  wt0 = which(BAF==0.0 & (!cnv.only))
  wt1 = which(BAF==1.0 & (!cnv.only))
  wt0NA = wt0[which(normalGtp[wt0] == -1)]
  wt0AA = wt0[which(normalGtp[wt0] ==  0)]
  wt0AB = wt0[which(normalGtp[wt0] ==  1)]
  wt1NA = wt1[which(normalGtp[wt1] == -1)]
  wt1BB = wt1[which(normalGtp[wt1] ==  2)]
  wt1AB = wt1[which(normalGtp[wt1] ==  1)]
  
 	for(z in 1:N){ 	  
    Hz  = H[z]
 	  bNz = matrix(0, nrow=L, ncol=Hz)
    idx = idxStart[z]

    if((z %in% c(2,4,6,8)) && contamination){
      pn0 = pnorm((0-mu.b[z,2])/sd.b[z,2])
      pn1 = pnorm((1-mu.b[z,3])/sd.b[z,3], lower.tail=FALSE)

      if(length(wt0NA) > 0){
        r0 = 0.5*ws[wt0NA,idx]
        r0 = r0/(r0 + pn0*ws[wt0NA,idx+1])
        bNz[wt0NA,1] = r0*pPN[wt0NA,z]
        bNz[wt0NA,2] = (1-r0)*pPN[wt0NA,z]
      }

      if(length(wt0AA) > 0){
        bNz[wt0AA,1] = pPN[wt0AA,z]
      }

      if(length(wt0AB) > 0){
        bNz[wt0AB,2] = pPN[wt0AB,z]
      }

      if(length(wt1NA) > 0){
        r1 = 0.5*ws[wt1NA,idx+3]
        r1 = r1/(r1 + pn1*ws[wt1NA,idx+2])
        bNz[wt1NA,4] = r1*pPN[wt1NA,z]
        bNz[wt1NA,3] = (1-r1)*pPN[wt1NA,z]
      }

      if(length(wt1BB) > 0){
        bNz[wt1BB,4] = pPN[wt1BB,z]
      }

      if(length(wt1AB) > 0){
        bNz[wt1AB,3] = pPN[wt1AB,z]
      }
    }else{
      bNz[wt0,1]  = pPN[wt0,z]
 	    bNz[wt1,Hz] = pPN[wt1,z]
    }

    for(hh in 1:Hz){
      ## for state 2/4/6/8, the 2nd/3rd components are not used
      ## if contamination is FALSE. 
      ## In that case the corresponding weights are 0
      ## so it does not hurt to do two extra summation
      ## but need to make sure the sds are all positive
      tmp  = (1-pi.b[z])*ws[Omiga,idx+hh-1]
      tmp  = tmp*dnorm(BAF[Omiga],mu.b[z,hh],sd.b[z,hh])
      bNz[Omiga,hh] = tmp*pPN[Omiga,z]/eBAF[Omiga,z]
    }

    x = idx
    y = x + Hz - 1
    postP36N[,x:y] = bNz
  }
  
  ## now genearate results for all the SNPs
  ## including those with missing values

  pP = matrix(nrow=L0, ncol=N)
  postP36 = matrix(nrow=L0, ncol=MM)
  
  if(length(wNA)>0){
    pP[-wNA,] = pPN
    postP36[-wNA,] = postP36N
  }else{
    pP = pPN
    postP36 = postP36N
  }
  
  # ------------------------------------------------
  ## posteriors of HMM states
  # ------------------------------------------------

  state.max = rep(NA, L0)
  if(length(wNA)>0){
    state.max[-wNA] = apply(pP[-wNA,], 1, which.max)
  }else{
    state.max = apply(pP, 1, which.max)
  }
  
  state.maxpP = apply(pP, 1, max)
  state.df    = data.frame(state=state.max, maxpP=state.maxpP)

  # ------------------------------------------------
  ## posteriors of copy number
  # ------------------------------------------------
  
  CN.pP = matrix(0, L0, 5)
  CN.pP[,1] = pP[,3]    # copy number: 0
  CN.pP[,2] = pP[,4]    # copy number: 1
  CN.pP[,3] = rowSums(pP[,1:2])  # copy number: 2

  if(CNA){
    CN.pP[,4] = rowSums(pP[,5:6])  # copy number: 3
    CN.pP[,5] = rowSums(pP[,7:9])  # copy number: 4
  }else{
    CN.pP[,4] = pP[,5]  # copy number: 3
    CN.pP[,5] = pP[,6]  # copy number: 4
  }
  colnames(CN.pP) = c("0", "1", "2", "3", "4")

  ## most likely copy number
  CN.max       = rep(NA, L0)
  if(length(wNA)>0){
    CN.max[-wNA] = apply(CN.pP[-wNA,],1,which.max) - 1
  }else{
    CN.max = apply(CN.pP,1,which.max) - 1
  }
  ## largest copy number post probability
  CN.maxpP = apply(CN.pP,1,max)
  CN.df    = data.frame(CN=CN.max, maxpP=CN.maxpP)

  # ------------------------------------------------
  ## posteriors of genotype classes
  # ------------------------------------------------

  Gtp.pP   = matrix(0, L0, 15)

  if(CNA){
    Gtp.pP[,1] = rowSums(postP36[,c(1,4,5)])          # "AA"
    Gtp.pP[,2] = postP36[,2]                          # "AB"
    Gtp.pP[,3] = rowSums(postP36[,c(3,6,7)])          # "BB"
    Gtp.pP[,4] = postP36[,8]                          # "Null"
    Gtp.pP[,5] = postP36[,9] + postP36[,10]           # "A"
    Gtp.pP[,6] = postP36[,11] + postP36[,12]          # "B"
    Gtp.pP[,7] = rowSums(postP36[,c(13,22,23)])       # "AAA"
    Gtp.pP[,8] = postP36[,14]                         # "AAB"
    Gtp.pP[,9] = postP36[,15]                         # "ABB"
    Gtp.pP[,10]= rowSums(postP36[,c(16,24,25)])       # "BBB"
    Gtp.pP[,11]= rowSums(postP36[,c(26,29,30,33)])    # "AAAA"
    Gtp.pP[,12]= postP36[,34]                         # "AAAB"
    Gtp.pP[,13]= postP36[,27]                         # "AABB"
    Gtp.pP[,14]= postP36[,35]                         # "ABBB"
    Gtp.pP[,15]= rowSums(postP36[,c(28,31,32,36)])    # "BBBB"    
  }else{
    Gtp.pP[,1] = rowSums(postP36[,c(1,4)])            # "AA"
    Gtp.pP[,2] = postP36[,2]                          # "AB"
    Gtp.pP[,3] = rowSums(postP36[,c(3,7)])            # "BB"
    Gtp.pP[,4] = postP36[,8]                          # "Null"
    Gtp.pP[,5] = postP36[,9]                          # "A"
    Gtp.pP[,6] = postP36[,12]                         # "B"
    Gtp.pP[,7:15] = postP36[,13:21]                   # "AAA" -- "BBBB"
  }
  
  colnames(Gtp.pP) = c("AA", "AB", "BB", "Null", "A", "B",
  "AAA", "AAB", "ABB", "BBB", "AAAA", "AAAB", "AABB", "ABBB", "BBBB")

  ## the most likely genotype
  Gtp.max       = rep(NA, L0)
  if(length(wNA) > 0){
    Gtp.max[-wNA] = apply(Gtp.pP[-wNA,],1,which.max)
  }else{
    Gtp.max = apply(Gtp.pP,1,which.max)
  }
  
  ## largest genotype post probability
  Gtp.maxpP = apply(Gtp.pP,1,max)
  Gtp.char  = colnames(Gtp.pP)[Gtp.max]
  w2drop    = which(Gtp.maxpP < 0.01)
  Gtp.char[w2drop]  = NA
  Gtp.maxpP[w2drop] = NA
  genotype.df = data.frame(genotype=Gtp.char, maxpP=Gtp.maxpP)
  
  postPL = list(state=pP, state.max=state.df, component=postP36, cn=CN.pP, 
    cn.max=CN.df, genotype=Gtp.pP, geno.max=genotype.df)
  
  class(postPL) = "postP.genoCN"
  postPL
}

