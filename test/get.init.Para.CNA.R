init.Para = list()

## parameters for LRR of SNP probes

init.Para[["pi.r"]] = c(0.01, 0.0, 0.0, 0.01, 0.01, 0, 0, 0, 0)
init.Para[["mu.r"]] = c(0.0,  0.0,  -3.5, -0.62, 0.37, 0.37, 0.65, 0.65, 0.65)
init.Para[["sd.r"]] = c(0.20, 0.20, 1.33, 0.20,  0.20, 0.20, 0.20, 0.20, 0.20)

init.Para[["mu.r.upper"]] = c(0.03, 0.03, -3, -0.25, 0.42, 0.42, 0.7, 0.7, 0.7)
init.Para[["mu.r.lower"]] = c(-.03, -.03, -4, -0.80, 0.30, 0.30, .58, .58, .58)
init.Para[["sd.r.upper"]] = c(0.25, 0.25, 1.5, 0.25, 0.25, 0.25, .25, .25, .25)
init.Para[["sd.r.lower"]] = c(0.15, 0.15, 1.2, 0.10, 0.15, 0.15, .15, .15, .15)

## parameters for LRR of CN-only probes

init.Para[["pi.rcn"]] = c(0.01, 0.0, 0.0, 0.01, 0.01, 0, 0, 0, 0)
init.Para[["mu.rcn"]] = c(0.0,  0.0,  -3.5, -0.62,  0.37, 0.37, 0.65, 0.65, 0.65)
init.Para[["sd.rcn"]] = c(0.20, 0.20, 1.33,  0.20,  0.20, 0.20, 0.20, 0.20, 0.20)

init.Para[["mu.rcn.upper"]] = c(0.03, 0.03, -3, -0.25, 0.42, 0.42, 0.7, 0.7, 0.7)
init.Para[["mu.rcn.lower"]] = c(-.03, -.03, -4, -0.80, 0.30, 0.30, .58, .58, .58)
init.Para[["sd.rcn.upper"]] = c(0.25, 0.25, 1.5, 0.25, 0.25, 0.25, .25, .25, .25)
init.Para[["sd.rcn.lower"]] = c(0.15, 0.15, 1.2, 0.10, 0.15, 0.15, .15, .15, .15)

init.Para[["pi.b"]] = rep(0.01, 9)

## matrix of mean parameter for BAF
## each row of mu.b corresponds to one state
mu.b = matrix(NA, nrow=9, ncol=4)
mu.b[1,1:3] = c(0, 0.5, 1)
mu.b[2,1:4] = c(0, 0.10, 0.90, 1)
mu.b[3,1]   = 0.5
mu.b[4,1:4] = c(0, 0.17, 0.83, 1)
mu.b[5,1:4] = c(0, 0.34, 0.66, 1)
mu.b[6,1:4] = c(0, 0.07, 0.93, 1)
mu.b[7,1:3] = c(0, 0.5, 1)
mu.b[8,1:4] = c(0, 0.055, 0.945, 1)
mu.b[9,1:4] = c(0, 0.25, 0.75, 1)

init.Para[["mu.b"]] = mu.b

## number of normal mixtures for states
H = c(3, 4, 1, 4, 4, 4, 3, 4, 4)

mu.b.upper = mu.b + 0.03
mu.b.lower = mu.b - 0.03
mu.b.upper[,1] = mu.b.lower[,1] = 0.0
for(z in 1:9){
  mu.b.upper[z,H[z]] = mu.b.lower[z,H[z]] = 1.0
}

mu.b.upper[2,2:3] = c(0.25, 0.98)
mu.b.lower[2,2:3] = c(0.02, 0.75)
mu.b.upper[4,2:3] = c(0.33, 0.98)
mu.b.lower[4,2:3] = c(0.02, 0.67)
mu.b.upper[6,2:3] = c(0.20, 0.98)
mu.b.lower[6,2:3] = c(0.02, 0.80)
mu.b.upper[8,2:3] = c(0.17, 0.98)
mu.b.lower[8,2:3] = c(0.02, 0.83)

mu.b.upper[3,1] = mu.b.lower[3,1] = 0.5
mu.b.upper[5,2:3] = c(0.40, 0.69)
mu.b.lower[5,2:3] = c(0.31, 0.60)
mu.b.upper[9,2:3] = c(0.33, 0.78)
mu.b.lower[9,2:3] = c(0.22, 0.67)

init.Para[["mu.b.upper"]] = mu.b.upper
init.Para[["mu.b.lower"]] = mu.b.lower

## matrix of sd parameter for BAF
## each row of sd.b corresponds to one state
sd.b = sd.b.lower=sd.b.upper=matrix(NA, nrow=9, ncol=4)
sd.b[1,1:3]       = c(.015, 0.03, .015)
sd.b.lower[1,1:3] = c(.008, .015, .008)
sd.b.upper[1,1:3] = c(0.03, 0.05, 0.03)
sd.b[3,1]         = 0.13
sd.b.lower[3,1]   = 0.10
sd.b.upper[3,1]   = 0.15
for(k in c(2,4,5,6,8,9)){
  sd.b[k,1:4]       = c(.015, 0.03, 0.03, .015)
  sd.b.lower[k,1:4] = c(.008, .015, .015, .008)
  sd.b.upper[k,1:4] = c(0.03, 0.05, 0.05, 0.03)
}

sd.b[7,1:3] = c(.015, 0.03, .015)
sd.b.lower[7,1:3] = c(.008, 0.015, .008)
sd.b.upper[7,1:3] = c(0.03, 0.05, 0.03)
init.Para[["sd.b"]] = sd.b
init.Para[["sd.b.upper"]] = sd.b.upper
init.Para[["sd.b.lower"]] = sd.b.lower

trans.m = matrix(1/8, nrow=9, ncol=9)
diag(trans.m)=0.0
init.Para[["trans.m"]] = trans.m
trans.begin = matrix(c(0.5, .05, .05, 0.1, 0.1, .05, .05, .05, .05), nrow=1)
init.Para[["trans.begin"]] = trans.begin

init.Para.CNA = init.Para
save(init.Para.CNA, file = "init.Para.CNA.rda")
