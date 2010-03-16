init.Para = list()

## parameters for LRR of SNP probes

init.Para[["pi.r"]] = c(0.01, 0.01, 0.10, 0.01, 0.01, 0.01)
init.Para[["mu.r"]] = c(0.0, 0.0, -3.5, -0.67, 0.40, 0.68)
init.Para[["sd.r"]] = c(0.15, 0.15, 1.3, 0.28, 0.20, 0.19)

init.Para[["mu.r.upper"]] = c(0.02, 0.02, -3.3, -0.665, 0.420, 0.700)
init.Para[["mu.r.lower"]] = c(-.02, -.02, -3.7, -0.680, 0.395, 0.675)
init.Para[["sd.r.upper"]] = c(0.175, 0.175, 1.4, 0.285, 0.205, 0.195)
init.Para[["sd.r.lower"]] = c(0.145, 0.145, 1.2, 0.260, 0.180, 0.180)

## parameters for LRR of CN-only probes

init.Para[["pi.rcn"]] = c(0.01, 0.01, 0.10, 0.01, 0.01, 0.01)
init.Para[["mu.rcn"]] = c(0.0,  0.0,  -2.1, -0.58,  0.37, 0.65)
init.Para[["sd.rcn"]] = c(0.18, 0.18,  2.1,  0.35,  0.24, 0.30)

init.Para[["mu.rcn.upper"]] = c(0.02, 0.02, -1.9, -0.575, 0.380, 0.70)
init.Para[["mu.rcn.lower"]] = c(-.02, -.02, -2.4, -0.590, 0.365, 0.63)
init.Para[["sd.rcn.upper"]] = c(0.210, 0.210, 2.2, 0.39, 0.25, 0.35)
init.Para[["sd.rcn.lower"]] = c(0.175, 0.175, 1.9, 0.30, 0.20, 0.20)

init.Para[["pi.b"]] = c(0.01, 0.01, 0.5, 0.01, 0.01, 0.01)

## matrix of mean parameter for BAF
## each row of mu.b corresponds to one state
mu.b = matrix(NA, nrow=6, ncol=5)
mu.b[1,1:3] = c(0, 0.5, 1)
## mu.b[2,2:3] will not be used since the corresponding weights are 0
mu.b[2,1:4] = c(0, NA, NA, 1) 
mu.b[3,1]   = 0.5
## mu.b[4,2:3] will not be used since the corresponding weights are 0
mu.b[4,1:4] = c(0, NA, NA, 1)
mu.b[5,1:4] = c(0, 0.333, 0.667, 1)
mu.b[6,1:5] = c(0, 0.25, 0.5, 0.75, 1)

init.Para[["mu.b"]] = mu.b

## number of normal mixtures for states
H = c(3, 4, 1, 4, 4, 5)

mu.b.upper = mu.b + 0.015
mu.b.lower = mu.b - 0.015
mu.b.upper[,1] = mu.b.lower[,1] = 0.0
for(z in 1:6){
  mu.b.upper[z,H[z]] = mu.b.lower[z,H[z]] = 1.0
}
mu.b.upper[3,1] = mu.b.lower[3,1] = 0.5

init.Para[["mu.b.upper"]] = mu.b.upper
init.Para[["mu.b.lower"]] = mu.b.lower

## matrix of sd parameter for BAF
## each row of sd.b corresponds to one state
sd.b = sd.b.lower=sd.b.upper=matrix(NA, nrow=6, ncol=5)
sd.b[1,1:3]       = c(0.016, 0.035, 0.016)
sd.b.lower[1,1:3] = c(0.005, 0.030, 0.005)
sd.b.upper[1,1:3] = c(0.020, 0.040, 0.020)
sd.b[3,1]         = 0.15
sd.b.lower[3,1]   = 0.10
sd.b.upper[3,1]   = 0.20
for(k in c(2,4)){
  sd.b[k,1:4]       = c(0.016, NA, NA, 0.016)
  sd.b.lower[k,1:4] = c(0.005, NA, NA, 0.005)
  sd.b.upper[k,1:4] = c(0.020, NA, NA, 0.020)
}

sd.b[5,1:4]       = c(0.016, 0.042, 0.042, 0.016)
sd.b.lower[5,1:4] = c(0.005, 0.025, 0.025, 0.005)
sd.b.upper[5,1:4] = c(0.020, 0.043, 0.043, 0.020)

sd.b[6,1:5] = c(0.016, 0.042, 0.035, 0.042, 0.016)
sd.b.lower[6,1:5] = c(0.005, 0.025, 0.025, 0.025, 0.005)
sd.b.upper[6,1:5] = c(0.020, 0.043, 0.038, 0.043, 0.020)

init.Para[["sd.b"]] = sd.b
init.Para[["sd.b.upper"]] = sd.b.upper
init.Para[["sd.b.lower"]] = sd.b.lower

transM = matrix(0, nrow=6, ncol=6)
transM[1,] = c(0.00, 0.01, 0.09, 0.80, 0.09, 0.01)
transM[2,] = c(0.20, 0.00, 0.20, 0.20, 0.20, 0.20)
transM[3,] = c(0.96, 0.01, 0.00, 0.01, 0.01, 0.01)
transM[4,] = c(0.96, 0.01, 0.01, 0.00, 0.01, 0.01)
transM[5,] = c(0.96, 0.01, 0.01, 0.01, 0.00, 0.01)
transM[6,] = c(0.96, 0.01, 0.01, 0.01, 0.01, 0.00)
init.Para[["trans.m"]] = transM

t1 = c(0.995, 0.005*c(.01, .09, .8, .09, .01))
trans.begin = matrix(t1/sum(t1), nrow=1)
init.Para[["trans.begin"]] = trans.begin

init.Para.CNV = init.Para
save(init.Para.CNV, file = "init.Para.CNV.rda")

