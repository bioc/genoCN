library(genoCN)

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

Theta = genoCNV(snpNames, chr, pos, LRR, BAF, pBs, 
            sampleID, cnv.only=cnv.only, outputSeg = TRUE, 
            outputSNP = 1, outputTag = "simu1")

seg = read.table("simu1_segment.txt", header=TRUE)
seg

snp = read.table("simu1_snp.txt", header=TRUE)
dim(snp)
snp[1:2,]

plotCN(pos, LRR, BAF, chr2plot = "22", fileNames="simu1_segment.txt")
