\name{init.Para.CNA}
\docType{data}
\alias{init.Para.CNA}
\title{ Initial parameters for the HMM }
\description{
  a list of initial values for the parameters of genoCNA. 
}
\usage{
data(init.Para.CNA)
}
\format{
  The format is a list of 16 items
\itemize{
\item{pi.r}{ a vector of length N, where N is the number of states.
   pi.r[j] is the prior probability of the uniform component of
   log R ratio for state j }
  \item{mu.r}{ a vector of length N, where N is the number of states.
   mu.r[j] is mean value of the normal component of
   log R ratio for state j }
  \item{sd.r}{ a vector of length N, where N is the number of states.
   sd.r[j] is standard deviation of the normal component of
   log R ratio for state j }
  \item{mu.r.upper, mu.r.lower}{ two vectors of the same size of mu.r,
   indicating the upper/lower bound of mu.r }
  \item{sd.r.upper, sd.r.lower}{ two vectors of the same size of sd.r,
   indicating the upper/lower bound of sd.r }
  \item{pi.b}{ a vector of length N, where N is the number of states.
   pi.b[j] is the prior probability of the uniform component of
   B allele frequency for state j }
  \item{mu.b}{ a matrix of N*M, where N is the number of states,
    and M is the maximum number of components of each states.
    mu.b[i,j] indicates the mean value of the j-th component of
    the i-th state }
  \item{sd.b}{ a matrix of the same size of mu.b, specifying the
    standard deviations }
  \item{mu.b.upper, mu.b.lower}{ two matrices of the same size of mu.b,
   incating the upper/lower bound of mu.b }
  \item{sd.b.upper, sd.b.lower}{ two matrices of the same size of sd.b,
   indicating the upper/lower bound of sd.b }
  \item{trans.m}{ transition probability matrix of size N*N.
  The diagonal elements are not used. }
  \item{trans.begin}{ a matrix of size S*N, where S is the number
  of chromosomes, and N is the number of states.  trans.begin[s,] are
  the state probabilities for the fist probe of the s-th chromosome.
  By default, we assume there is only one chromosome, therefore it is
  a matrix of 1*N.}
}}
\examples{
data(init.Para.CNA)
}
\keyword{datasets}
