
# ----------------------------------------------------------------------
# check parameters
# ----------------------------------------------------------------------

check.para <- function(Para, CNA, nChr)
{
  # ------------------------------------------------
  # parameters to check
  # ------------------------------------------------

  mu.b = Para[["mu.b"]]
  sd.b = Para[["sd.b"]]
  mu.b.lower = Para[["mu.b.lower"]]
  mu.b.upper = Para[["mu.b.upper"]]
  sd.b.lower = Para[["sd.b.lower"]]
  sd.b.upper = Para[["sd.b.upper"]]
  
  sd.r = Para[["sd.r"]]
  trans.begin = Para[["trans.begin"]]
  
  # ------------------------------------------------
  # check parameters
  # ------------------------------------------------
  if(any(sd.r <= 0)){
    stop("sd.r must be greater than 0\n")
  }
  
  if(is.vector(trans.begin)){
  	N = length(trans.begin)
    trans.begin=matrix(trans.begin, nrow=nChr, ncol=N, byrow=TRUE)
  }else{
  	if(nrow(trans.begin) == 1){
  	  trans.begin = as.vector(trans.begin)
  	  N = length(trans.begin)
      trans.begin=matrix(trans.begin, nrow=nChr, ncol=N, byrow=TRUE)
    }
    if(nrow(trans.begin) != nChr){
      stop("number of rows of trans.begin does not match the number of chromosomes\n")
    }	
  }
  
  ## number of states
  N = ncol(trans.begin)
  
  if(!CNA && N!=6){
    stop("don't we have 6 states for CNV studies?")
  }

  if(CNA && N!=9){
    stop("don't we have 9 states for CNA studies?")
  }

  # ------------------------------------------------
  # check mu.b, mu.b.lower, and mu.b.upper
  # ------------------------------------------------
  if(!CNA){
    mu.b.use = matrix(FALSE, nrow=6, ncol=5)
    mu.b.use[1,1:3]    = TRUE
    mu.b.use[2,c(1,4)] = TRUE
    mu.b.use[3,1]      = TRUE
    mu.b.use[4,c(1,4)] = TRUE
    mu.b.use[5,1:4]    = TRUE
    mu.b.use[6,1:5]    = TRUE
  }else{
    mu.b.use = matrix(FALSE, nrow=9, ncol=4)
    mu.b.use[1,1:3] = TRUE
    mu.b.use[2,1:4] = TRUE
    mu.b.use[3,1]   = TRUE
    mu.b.use[4,1:4] = TRUE
    mu.b.use[5,1:4] = TRUE
    mu.b.use[6,1:4] = TRUE
    mu.b.use[7,1:3] = TRUE
    mu.b.use[8,1:4] = TRUE
    mu.b.use[9,1:4] = TRUE
  }
  
  if(any(dim(mu.b.upper) != dim(mu.b))){
    stop("dimension of mu.b.upper does not match with mu.b\n")
  }

  if(any(dim(mu.b.lower) != dim(mu.b))){
    stop("dimension of mu.b.lower does not match with mu.b\n")
  }

  mu.b.vect = as.numeric(mu.b[mu.b.use])
  if(any(is.na(mu.b.vect))){
    stop("some neccesary parameters in mu.b are missing")
  }

  mu.b.upper.vect = as.numeric(mu.b.upper[mu.b.use])
  if(any(is.na(mu.b.upper.vect))){
    stop("some neccesary parameters in mu.b.upper are missing")
  }

  mu.b.lower.vect = as.numeric(mu.b.lower[mu.b.use])
  if(any(is.na(mu.b.lower.vect))){
    stop("some neccesary parameters in mu.b.lower are missing")
  }

  # ------------------------------------------------
  # check sd.b, sd.b.lower, and sd.b.upper
  # ------------------------------------------------

  if(any(dim(sd.b.upper) != dim(sd.b))){
    stop("dimension of sd.b.upper does not match with mu.b\n")
  }

  if(any(dim(sd.b.lower) != dim(sd.b))){
    stop("dimension of sd.b.lower does not match with mu.b\n")
  }

  sd.b.use  = mu.b.use
  sd.b.vect = as.numeric(sd.b[sd.b.use])
  if(any(is.na(sd.b.vect))){
    stop("some neccesary parameters in sd.b are missing\n")
  }
  if(any(sd.b.vect <= 0)){
    stop("sd.b must be greater than 0\n")
  }
  
  sd.b.upper.vect = as.numeric(sd.b.upper[sd.b.use])
  if(any(is.na(sd.b.upper.vect))){
    stop("some neccesary parameters in sd.b.upper are missing\n")
  }
  if(any(sd.b.upper.vect <= 0)){
    stop("sd.b.upper must be greater than 0\n")
  }

  sd.b.lower.vect = as.numeric(sd.b.lower[sd.b.use])
  if(any(is.na(sd.b.lower.vect))){
    stop("some neccesary parameters in sd.b.lower are missing\n")
  }
  if(any(sd.b.lower.vect <= 0)){
    stop("sd.b.lower must be greater than 0\n")
  }

  mu.b[is.na(mu.b)] = 0.0
  sd.b[is.na(sd.b)] = 0.0
  mu.b.lower[is.na(mu.b.lower)] = 0.0
  mu.b.upper[is.na(mu.b.upper)] = 0.0
  sd.b.lower[is.na(sd.b.lower)] = 0.0
  sd.b.upper[is.na(sd.b.upper)] = 0.0

  Para[["mu.b"]] = mu.b
  Para[["sd.b"]] = sd.b
  Para[["mu.b.lower"]] = mu.b.lower
  Para[["mu.b.upper"]] = mu.b.upper
  Para[["sd.b.lower"]] = sd.b.lower
  Para[["sd.b.upper"]] = sd.b.upper
  
  Para[["sd.r"]] = sd.r
  Para[["trans.begin"]] = trans.begin
  
  Para
}
