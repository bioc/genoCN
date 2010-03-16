\name{snpData}
\alias{snpData}
\docType{data}
\title{
  Simulated LRR and BAF data for 17,348 SNPs on chromosome 22. 
}
\description{
   Simulated LRR and BAF data for 17,348 SNPs on chromosome 22. Two CNVs are simulated. One is from the 1001-th probe to the 1100-th probe, with copy number 1. The other one is from the 10,001-th probe to the 10,200-th probe, with copy number 3. 
}
\usage{data(snpData)}
\format{
  A data frame with 17,348 observations on the following 3 variables.
  \describe{
    \item{\code{Name}}{a character vector of probe Names}
    \item{\code{LRR}}{a numeric vector of LRR values of each probe}
    \item{\code{BAF}}{a numeric vector of BAF of each probe}
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
