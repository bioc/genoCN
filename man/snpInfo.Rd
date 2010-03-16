\name{snpInfo}
\alias{snpInfo}
\docType{data}
\title{
  Information of 17,348 SNPs on chromosome 22.
}
\description{
  Information of 17,348 SNPs on chromosome 22.
}
\usage{data(snpInfo)}
\format{
  A data frame with 17348 observations on the following 4 variables.
  \describe{
    \item{\code{Name}}{a character vector of probe Names}
    \item{\code{Chr}}{a character vector of chromosomes of each probe}
    \item{\code{Position}}{a numeric vector of genomic position of each probe}
    \item{\code{PFB}}{a numeric vector of population frequency of B allele for each probe. For copy number only probes, PFB=2.0}
  }
}
\examples{

data(snpData)
data(snpInfo)

dim(snpData)
dim(snpInfo)

snpData[1:2,]
snpInfo[1:2,]

plotCN(pos=snpInfo$Position, LRR=snpData$LRR, BAF=snpData$BAF, 
main = "simulated data on Chr22")

}
\keyword{datasets}
