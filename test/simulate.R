library(genoCN)
setwd("~/research/cfcnv/UNC/")

info = read.table("snp_info_hg18.txt", header=TRUE, sep="\t", as.is=TRUE)
dim(info)
info[1:2,]

snpInfo = info[info$Chr=="22", ]
dim(snpInfo)
snpInfo[1:2,]

#--------------------------------------------------------------------
# simulate the underlying states:
#  100 probes with copy number 1
#  200 probes with copy number 3
#--------------------------------------------------------------------

n = nrow(snpInfo)
z = rep(1, n)

z[1001:1100] = 4
z[10001:10200] = 5

#--------------------------------------------------------------------
# simulate LRR
#--------------------------------------------------------------------

r = rnorm(n, 0, 0.15)
r[1001:1100] = rnorm(100, -0.67, 0.28)
r[10001:10200] = rnorm(200, 0.40, 0.20)

#--------------------------------------------------------------------
# simulate genotype and BAF when copy number is 2
#--------------------------------------------------------------------

geno = b = rep(NA, n)
wuse = which(snpInfo$PFB <=1)

for(w in wuse){
  geno[w] = rbinom(1, 2, snpInfo$PFB[w])
}
table(geno)

w0 = which(geno==0)
w1 = which(geno==1)
w2 = which(geno==2)

bw0 = rnorm(length(w0), 0.0, 0.016)
bw1 = rnorm(length(w1), 0.5, 0.035)
bw2 = rnorm(length(w2), 1.0, 0.016)

bw0[which(bw0 < 0)] = 0
bw2[which(bw2 > 1)] = 1
bw1[which(bw1 < 0)] = 0
bw1[which(bw1 > 1)] = 1

b[w0] = bw0
b[w1] = bw1
b[w2] = bw2

#--------------------------------------------------------------------
# simulate genotype and BAF when copy number is 1
#--------------------------------------------------------------------

geno1 = rep(NA, 100)
for(i in 1:100){
  geno1[i] = rbinom(1, 1, snpInfo$PFB[1000+i])
}
table(geno1)

w0 = which(geno1==0)
w1 = which(geno1==1)

bw0 = rnorm(length(w0), 0.0, 0.016)
bw1 = rnorm(length(w1), 1.0, 0.016)

bw0[which(bw0 < 0)] = 0
bw1[which(bw1 > 1)] = 1

b[w0+1000] = bw0
b[w1+1000] = bw1


#--------------------------------------------------------------------
# simulate genotype and BAF when copy number is 3
#--------------------------------------------------------------------

geno3 = rep(NA, 200)
for(i in 1:200){
  geno3[i] = rbinom(1, 3, snpInfo$PFB[10000+i])
}
table(geno3)

w0 = which(geno3==0)
w1 = which(geno3==1)
w2 = which(geno3==2)
w3 = which(geno3==3)

bw0 = rnorm(length(w0), 0.0,   0.016)
bw1 = rnorm(length(w1), 0.333, 0.042)
bw2 = rnorm(length(w2), 0.667, 0.042)
bw3 = rnorm(length(w3), 1.0,   0.016)

bw0[which(bw0 < 0)] = 0
bw3[which(bw3 > 1)] = 1

b[w0+10000] = bw0
b[w1+10000] = bw1
b[w2+10000] = bw2
b[w3+10000] = bw3


hist(b, breaks=100)

hist(r, breaks=100)

plotCN(fileNames=NULL, snpInfo$Position, LRR=r, BAF=b)

snpData = data.frame(Name=snpInfo$Name, LRR=r, BAF=b)
snpData$Name = as.character(snpData$Name)

save(snpData, file="snpData.rda")
save(snpInfo, file="snpInfo.rda")



