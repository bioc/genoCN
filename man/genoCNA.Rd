\name{genoCNA}
\alias{genoCNA}
\title{ Copy Number Aberration  }
\description{
  extract genotype and copy number calls for copy number aberrations, which are often observed in tumor tissues
}
\usage{
genoCNA(snpNames, chr, pos, LRR, BAF, pBs, sampleID, 
  Para=NULL, fixPara=FALSE, cnv.only=NULL, estimate.pi.r=TRUE, 
  estimate.pi.b=TRUE, estimate.trans.m=TRUE, outputSeg = TRUE, 
  outputSNP=3, outputTag=sampleID, outputViterbi=FALSE, 
  Ds=c(1e10, 1e10, rep(1e8, 7)), pBs.alpha=0.001, contamination=TRUE, 
  normalGtp=NULL, geno.error=0.01, min.tp=1e-4, max.diff=0.1, 
  distThreshold=1e6, transB=c(0.5,.05,.05,0.1,0.1,.05,.05,.05,.05), 
  epsilon=0.005, K=5, maxIt=200, seg.nSNP=3, traceIt=5)
}
\arguments{
  \item{snpNames}{ a vector of SNP names. SNPs must be ordered by chromosme locations }
  
  \item{chr}{ chromosomes of all the SNPs specified in \code{snpNames} }
  
  \item{pos}{ positions of all the SNPs specified in \code{snpNames}}
  
  \item{LRR}{ Log R Ratio of all the SNPs specified in \code{snpNames}}
  
  \item{BAF}{ B Allele Frequency of all the SNPs specified in \code{snpNames}}
  
  \item{pBs}{ population frequency of of all the SNPs specified in \code{snpNames} }
  
  \item{sampleID}{ symbol/name of the studied sample. Only one sample is studied each time }
  
  \item{Para}{ a list of initial parameters for the HMM. If Para is NULL, The 
   default initial parameters: init.Para.CNA is used }

  \item{fixPara}{ if fixPara is TRUE, the parameters in Para are fixed, and are 
  used directly to calculate posterior probabilities. It is not recommended to 
  set fixPara as TRUE for CNA studies. }

  \item{cnv.only}{ a vector indicating those CNV-only probes, for which we 
   only consider their Log R ratio. If it is NULL, there is no CNV-only probes }

  \item{estimate.pi.r}{ to estimate pi.r (proportion of uniform component for 
  LRR) or not. By default, estimate.pi.r=FALSE, and the initial value of pi.r is 
  used to estimate other parameters }

  \item{estimate.pi.b}{ to estimate pi.b (proportion of uniform component for 
  BAF) or not. By default, estimate.pi.b=FALSE, and the initial value of pi.b is 
  used to estimate other parameters }
	
  \item{estimate.trans.m}{ to estimate transition probability matrix or not. By 
  default, estimate.trans.m=FALSE, and the initial value of estimate.trans.m is 
  used to estimate other parameters }
	
  \item{outputSeg}{ wether to output the information of copy number altered 
  segments }
  
  \item{outputSNP}{ if \code{outputSNP} is 0, do not output SNP specific information; if \code{outputSNP} is 1, output the most likely copy number and genotype state of the SNPs that are within copy number altered regions; if \code{outputSNP} is 2, output the most likely copy number and genotype state of all the SNPs (whether it is within CNV regions or not), if \code{outputSNP} is 3, output the posterior probability for all the copy number and genotype states for the SNPs. }
  
  \item{outputTag}{ the prefix of the output files, output of copy number 
  altered segments is written into file outputTag\_segment.txt, and output of 
  SNP information is written into file outputTag\_SNP.txt }
  
  \item{outputViterbi}{ whether to output the copy altered regions identified by 
  the viterbi algorithm. see details }
  
  \item{Ds}{ Parameter to for transition probability of the HMM. A vector of 
  length N, where N is the number of states in the HMM }
  
  \item{pBs.alpha}{
    pBs.alpha is the lower limit of population B allele frequency, and the 
    upper limit is 1 - pBs.alpha
  }
  
  \item{contamination}{ whether tissue contamination is considered }
  
  \item{normalGtp}{ \code{normalGtp} is specified only if paired tumor-normal 
  SNP array is availalble. It is the normal tissue genotype for all the SNPs 
  specified in \code{snpNames}, which can only take four different values:
  -1, 0, 1, and 2. Values 0, 1, 2 correspond to the number of B alleles, 
  and value -1 indicates the normal genotype is missing. By default, 
  it is NULL, then all the normal genotype are set missing (-1) }
  
  \item{geno.error}{ probability of genotyping error in normal tissue genotypes }
  
  \item{min.tp}{ the minimum of transition probability. }
  
  \item{max.diff}{ Due to normalization procedure, the BAF may not be symmetric. 
  Let's use state (AAA, AAB, ABB, BBB) as an example. Ideally, mean values of 
  normal components AAB and ABB, denoted by mu1 and mu2, respectively, should 
  have the relation mu1 = 1-mu2 if BAF is symmetric. However, this may not be 
  true due to normalization procedures. We restrict the difference of mu1 and 
  (1-mu2) by this parameter max.diff. }

  \item{distThreshold}{ If distance between adjacent probes is larger than 
  distThreshold, restart the transition probability by the default values
  in \code{transB}.}

  \item{transB}{ The default transition probability. }
  
  \item{epsilon}{ see explanation of \code{K}}
  
  \item{K}{ epsilon and K are used to specify the convergence
  criteria. We say the estimate.para is converged if for K consecutive 
  updates, the maximum change of parameter estimates in every adjacent
  step is smaller than epsilon }
  
  \item{maxIt}{
    the maximum number of iterations of the EM algorithm to estimate parameters
  }

  \item{seg.nSNP}{
    the minimum number of SNPs per segment
  }
  
  \item{traceIt}{
    if traceIt is a integer n, then the running time is printed out in every n iterations of the EM algorithm.
    if traceIt is 0 or negative, no tracing information is printed out.
  }
}
\value{
  results are written into output files
}
\note{
  Copy number altered regions are identified, by default, based on the SNP level 
  copy number calls. A CNA region boundary is declared simply when the adjacent 
  SNPs have different copy numbers. An alternative approach is to use viterbi 
  algorithm to output the ``best path''. Most time the results based on the SNP 
  level copy number calls are the same as the results from viterbi algorithm. 
  For the following up association studies, the SNP level information is more 
  relevant if we examine the association SNP by SNP. 
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

snpNames = snpInfo$Name
chr = snpInfo$Chr
pos = snpInfo$Position
LRR = snpData$LRR
BAF = snpData$BAF
pBs = snpInfo$PFB
cnv.only=(snpInfo$PFB>1)
sampleID="simu1"

# Note this simulated data is more of CNV rather than CNA. 
# For example, there is no tissue contamination. 
# We just use it to illustrate the usage of genoCNA. 

Theta = genoCNA(snpNames, chr, pos, LRR, BAF, pBs, contamination=TRUE, 
  normalGtp=NULL, sampleID, cnv.only=cnv.only, outputSeg = TRUE, 
            outputSNP = 1, outputTag = "simu1")
}
\author{ Wei Sun and Zhengzheng Tang }
\keyword{ methods }
