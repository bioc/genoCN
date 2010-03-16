# ----------------------------------------------------------------------
#  log sum exp
# ----------------------------------------------------------------------
logsumexp <- function(v){
  if(any(v==Inf)){
    stop("positive inifinite value in v\n")
  }else{
    w = which(v == -Inf)
    if(length(w)>0){
      v = v[-w]
    }
  }
  if(length(v)==0){
    lse = -Inf
  }else if(length(v)==1){
    lse = v[1]
  }else{
    wv  = which.max(v)
    mv  = max(v)
    res = sum(exp(v[-wv] - mv))
    lse = mv + log(1+res)
  }
  lse
}

# ----------------------------------------------------------------------
#  Maximum Likelihood Estimation of truncated normal distribution.
#  given the normal distribution is truncated in one side, and
#  we know the number of truncated observations.
# ----------------------------------------------------------------------

tnorm.mle <- function(x, Tr, n, roundoff=1e-6, maxIt=100){

  if(all(x < Tr)){
    left = 0
  }else if(all(x > Tr)){
    left = 1
  }else{
    stop("x must be all smaller than Tr or be all bigger than Tr\n")
  }

  M1 = mean(x)
  M2 = sum(x*x)/length(x)
  
  dims = numeric(3)
  dims[1] = length(x)
  dims[2] = n
  dims[3] = maxIt
  
  mu = sigma = 0;

  Z = .C("tnorm_mleR", as.double(Tr), as.double(M1), as.double(M2),
    as.integer(dims), mu=as.double(mu), sigma=as.double(sigma),
    as.integer(left), as.double(roundoff), PACKAGE="genoCN")
  list(mu = Z[["mu"]], sigma= Z[["sigma"]])
}

