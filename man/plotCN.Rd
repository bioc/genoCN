\name{plotCN}
\alias{plotCN}
\title{ plot LRR, BAF, and the copy number estimates }
\description{
   plot LRR, BAF, and the copy number estimates of genoCNV and/or PennCNV.
}
\usage{
plotCN(pos, LRR, BAF, chr2plot = NULL, sampleIDs = NULL, fileNames=NULL, 
types = "genoCN", CNA = TRUE, main = "", LRR.ylim=NULL,
cex=0.5, plot.lowess=TRUE)
}
\arguments{
  \item{pos}{ position of all the SNPs }
  \item{LRR}{ a vector of the log R ratio, should be one-to-one
  correspondence of pos }
  \item{BAF}{ a vector of the B allele frequency, should be one-to-one
  correspondence of pos }
  \item{chr2plot}{ which chromosome to plot. Only one chromosome can be plotted
  each time }
  \item{sampleIDs}{ sample ID, could be a vector of the same length as fileNames so that different sample IDs are used for different input files. }
  \item{fileNames}{ one or more names of the output files of genoCN or PennCNV. If it is NULL, only plot the LRR and BAF. }
  \item{types}{ should be the same length as fileNames, indicating the type of output, currently only support "genoCN" and "pennCNV" }
  \item{CNA}{ whether this is a copy number aberration study. }
  \item{main}{ title of the plot }
  \item{LRR.ylim}{ Range of y-axis for LRR plot }
  \item{cex}{ the amount by which plotting text
          and symbols should be magnified relative to the default}
  \item{plot.lowess}{ to plot the lowess curve for LRR or not }
}
\examples{
data(snpData)
data(snpInfo)

dim(snpData)
dim(snpInfo)

snpData[1:2,]
snpInfo[1:2,]

snpInfo[c(1001,1100,10001,10200),]

plotCN(pos=snpInfo$Position, LRR=snpData$LRR, BAF=snpData$BAF, 
main = "simulated data on Chr22")
}
\author{ Wei Sun }
\seealso{ \code{\link{genoCNA}}, \code{\link{genoCNV}} }
\keyword{ methods }
