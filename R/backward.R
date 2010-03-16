# ----------------------------------------------------------------------
#  return a backward probablility matrix in log scale
# ----------------------------------------------------------------------

backward <- function(pos, em, trans.m, Ds, distThreshold, transB){

  ## numer of states
  N = nrow(trans.m)

  ## sequence length
  L = length(pos)
      
  ## backward probability
  b = matrix(0, nrow=L, ncol=N)
  
  dims = c(L, N)
  
  em.inf = which(t(em)==-Inf)
  em.noinf = length(em.inf)  
  em[em==-Inf] = 0
  
  Z = .C("backward", as.double(pos), as.double(t(em)), as.integer(em.inf), 
         as.integer(em.noinf),as.double(t(trans.m)),  
         as.double(Ds), b=as.double(t(b)), as.integer(dims), 
         as.double(distThreshold), as.double(transB), PACKAGE="genoCN")
  
  b = matrix(Z$b, byrow=TRUE, nrow=L, ncol=N)
  
  b
}

# ----------------------------------------------------------------------
#  return a forward probablility matrix in log scale
# ----------------------------------------------------------------------

forward <- function(pos, em, trans.m, trans.begin, Ds, distThreshold, transB){
  
  ## numer of states
  N = nrow(trans.m)
  
  ## sequence length
  L = length(pos)
  
  ## forward probability matrix
  f = matrix(0, nrow=L, ncol=N)
  
  dims = c(L, N)
  
  em.inf = which(t(em)==-Inf)
  em.noinf = length(em.inf)  
  em[em==-Inf] = 0
  
  Z = .C("forward", as.double(pos), as.double(t(em)), as.integer(em.inf), 
         as.integer(em.noinf), as.double(t(trans.m)), as.double(trans.begin), 
         as.double(Ds), f=as.double(t(f)), as.integer(dims), 
         as.double(distThreshold), as.double(transB), PACKAGE="genoCN")
  f = matrix(Z$f, byrow=TRUE, nrow=L, ncol=N)
  
  f
}

# ----------------------------------------------------------------------
# posterior probability within forward-backward algorithm
# ----------------------------------------------------------------------

postP <- function(f, b){
  L = nrow(f)
  N = ncol(f)
  dims=c(L,N)
  if(nrow(b)!=L || ncol(b)!=N){
    stop("dimension of f and b do not match\n")
  }
  
  pP = matrix(0, nrow=L, ncol=N)
  
  f.inf = which(as.vector(t(f))==-Inf)
  f.noinf = length(f.inf)  
  f[f==-Inf] = 0
  b.inf = which(t(b)==-Inf)
  b.noinf = length(b.inf)
  b[b==-Inf] = 0
  
  no.inf = c(f.noinf, b.noinf)
  
  Z = .C( "postP", as.double(t(f)), as.double(t(b)), as.integer(f.inf), 
         as.integer(b.inf), as.integer(no.inf), as.integer(dims),
         pP=as.double(t(pP)), PACKAGE="genoCN" )
  
  pP = matrix(Z$pP, byrow=TRUE, nrow=L, ncol=N)
  
  xx = apply(pP, 1, sum)
  if(min(xx, na.rm=TRUE) < 0.999 || max(xx, na.rm=TRUE) > 1.001){
    stop("sum of posterior P is not 1\n")
  }
  
  pP
}

