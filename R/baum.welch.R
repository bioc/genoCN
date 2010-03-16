# ----------------------------------------------------------------------
#  one step update of baum.welch algorithm
# ----------------------------------------------------------------------

baum.welch <- function(pos, LRR, BAF, pBs, Para, R.m, Ds, cnv.only, 
  estimate.pi.r, estimate.pi.b, estimate.trans.m, CNA, min.tp, 
  max.diff, len, distThreshold, contamination, normalGtp, geno.error){

  ##
  if((!CNA) && contamination){
    stop("contamination is TRUE for CNV studies\n")
  }

  trans.begin = Para[["trans.begin"]]
  transB = Para[["transB"]]

  pi.r = Para[["pi.r"]]
  mu.r = Para[["mu.r"]]
  sd.r = Para[["sd.r"]]
  mu.r.upper = Para[["mu.r.upper"]]
  mu.r.lower = Para[["mu.r.lower"]]
  sd.r.upper = Para[["sd.r.upper"]]
  sd.r.lower = Para[["sd.r.lower"]]
  
  pi.rcn = Para[["pi.rcn"]]
  mu.rcn = Para[["mu.rcn"]]
  sd.rcn = Para[["sd.rcn"]]
  mu.rcn.upper = Para[["mu.rcn.upper"]]
  mu.rcn.lower = Para[["mu.rcn.lower"]]
  sd.rcn.upper = Para[["sd.rcn.upper"]]
  sd.rcn.lower = Para[["sd.rcn.lower"]]

  pi.b = Para[["pi.b"]]
  mu.b = Para[["mu.b"]]
  sd.b = Para[["sd.b"]]
  trans.m = Para[["trans.m"]]
  mu.b.upper = Para[["mu.b.upper"]]
  mu.b.lower = Para[["mu.b.lower"]]
  sd.b.upper = Para[["sd.b.upper"]]
  sd.b.lower = Para[["sd.b.lower"]]
  
  # ------------------------------------------------
  # dimensions
  # ------------------------------------------------
  ## sequence length
  L = length(pos)

  ## number of states
  N = ncol(trans.begin)

  ## maximum number of normal components in mixture distribution
  M = ncol(mu.b)

  dims = c(L, N, M)

  # ------------------------------------------------
  # emission probability
  # ------------------------------------------------
  
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
  # forward/backward probability
  # ------------------------------------------------
  S = length(len)
  
  x = 1
  f = b = pP = matrix(0, nrow=L, ncol=N)
  trans.begin.chr = numeric(N)
  logL = numeric(S) # overall log liklihood
    
  for(i in 1:S){
    f.chr  = matrix(0, nrow=len[i], ncol=N)
    b.chr  = matrix(0, nrow=len[i], ncol=N)
    pP.chr = matrix(0, nrow=len[i], ncol=N)
    
    y = x + len[i] - 1
    
    trans.begin.chr = trans.begin[i,]

    ## forward probability matrix
    f.chr = forward(pos[x:y], em[x:y,], trans.m, trans.begin.chr, 
                    Ds, distThreshold, transB)
  
    ## bakward probability matrix
    b.chr = backward(pos[x:y], em[x:y,], trans.m, Ds, 
                     distThreshold, transB)
  
    ## posterior probability matrix
    pP.chr = postP(f.chr, b.chr)
    
    f[x:y,]  = f.chr
    b[x:y,]  = b.chr
    pP[x:y,] = pP.chr
    
    x = x + len[i]
    logL[i] = logsumexp(f[y,])
  }

  logL = logsumexp(logL)
  
  rm(f.chr, b.chr, pP.chr)
  
  # ------------------------------------------------
  # one iteration of Baum-Welch algorithm
  # ------------------------------------------------
    
  pi.r.new = pi.r
  mu.r.new = numeric(N)
  sd.r.new = numeric(N)

  pi.rcn.new = pi.rcn
  mu.rcn.new = numeric(N)
  sd.rcn.new = numeric(N)

  pi.b.new = pi.b
  
  mu.b.new = matrix(0, nrow=N, ncol=M)
  sd.b.new = matrix(0, nrow=N, ncol=M)
    
  trans.m.new = trans.m
  trans.begin.new = trans.begin

  em.inf   = which(t(em)==-Inf)
  em.noinf = length(em.inf)
  em[em==-Inf] = 0

  f.inf   = which(t(f)==-Inf)
  f.noinf = length(f.inf)
  f[f==-Inf] = 0
  
  b.inf   = which(t(b)==-Inf)
  b.noinf = length(b.inf)
  b[b==-Inf] = 0
  
  no.inf = c(em.noinf, f.noinf, b.noinf)

  DEBUG = 0
  dims  = numeric(7)
  dims[1:6] = c(L, N, M, S, DEBUG, as.integer(estimate.pi.r))
  dims[7:8] = as.integer(c(estimate.pi.b, estimate.trans.m))

  Z = .C("baum_welch", as.double(pos), as.double(LRR), as.double(BAF),
    as.double(pBs), as.integer(normalGtp), as.double(geno.error), 
    as.integer(!cnv.only), as.integer(CNA), as.double(Ds), 
    as.double(R.m), as.double(pi.r), as.double(pi.rcn), 
    as.double(mu.r), as.double(mu.r.upper), as.double(mu.r.lower),
    as.double(sd.r), as.double(sd.r.upper), as.double(sd.r.lower),
    as.double(mu.rcn), as.double(mu.rcn.upper), as.double(mu.rcn.lower),
    as.double(sd.rcn), as.double(sd.rcn.upper), as.double(sd.rcn.lower),
    as.double(pi.b), as.double(t(trans.m)), as.double(min.tp),
    as.double(max.diff), as.double(t(mu.b)), as.double(t(mu.b.upper)),
    as.double(t(mu.b.lower)), as.double(t(sd.b)),
    as.double(t(sd.b.upper)), as.double(t(sd.b.lower)),
    as.double(t(em)), as.integer(em.inf), as.double(t(eLRR)),
    as.double(t(eBAF)), as.double(t(f)), as.double(t(b)),
    as.integer(f.inf), as.integer(b.inf), as.integer(no.inf),
    as.double(t(pP)), as.integer(contamination), 
    as.integer(dims), as.integer(len), as.double(distThreshold), 
    pi.r.new=as.double(pi.r.new), mu.r.new=as.double(mu.r.new),
    sd.r.new=as.double(sd.r.new), pi.rcn.new=as.double(pi.rcn.new), 
    mu.rcn.new=as.double(mu.rcn.new), sd.rcn.new=as.double(sd.rcn.new), 
    pi.b.new=as.double(pi.b.new),
    mu.b.new=as.double(t(mu.b.new)), sd.b.new=as.double(t(sd.b.new)),
    trans.m.new=as.double(t(trans.m.new)), 
    trans.begin.new=as.double(t(trans.begin.new)),
    PACKAGE="genoCN")
    
  mu.b.new = matrix(Z$mu.b.new, byrow=TRUE, nrow=N, ncol=M)
  sd.b.new = matrix(Z$sd.b.new, byrow=TRUE, nrow=N, ncol=M)
  trans.m.new = matrix(Z$trans.m.new, byrow=TRUE, nrow=N, ncol=N)
  trans.begin.new = matrix(Z$trans.begin.new, byrow=TRUE, nrow=S, ncol=N)
  
  Theta1 = list(pi.r=Z$pi.r.new, mu.r=Z$mu.r.new, sd.r=Z$sd.r.new,
                pi.rcn=Z$pi.rcn.new, mu.rcn=Z$mu.rcn.new, sd.rcn=Z$sd.rcn.new,
                pi.b=Z$pi.b.new, mu.b=mu.b.new, sd.b=sd.b.new,
                trans.m=trans.m.new, trans.begin=trans.begin.new,
                transB=transB, logL=logL)
  rm(Z)
  
  Theta1
}
