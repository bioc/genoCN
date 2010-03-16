`tnorm.mle` <-
function(x, Tr, n){

  if(all(x < Tr)){
    truncation = "right"
  }else if(all(x > Tr)){
    truncation = "left"
  }else{
    stop("x must be all smaller than Tr or be all bigger than Tr\n")
  }

  r = length(x)
  M1 = mean(x)
  M2 = sum(x^2)/r
  
  Vpn2 = 4*(M2 - 2*Tr*M1 + Tr*Tr)/(Tr - M1)/(Tr - M1)
  p = r/n

  gh = function(hh, p, Vpn2){
    gg = (p/((1-p)*Vpn2))*((-Vpn2+2)*hh + 2*sqrt(hh^2 + Vpn2))
    gg = gg - dnorm(hh)/pnorm(hh,lower.tail=FALSE)
    gg
  }

  ur = uniroot(gh, p=p, Vpn2=Vpn2, lower=-5, upper=5)
  hh = ur$root
  
  if(truncation == "left"){
    hh = -hh
    sigma = 0.5*(Tr - mean(x))*(-hh - sqrt(hh^2 + Vpn2))
  }else{
    sigma = 0.5*(Tr - mean(x))*(-hh + sqrt(hh^2 + Vpn2))
  }
  mu = Tr - sigma*hh
  list(h=hh, mu=mu, sigma=sigma)
}

n = 1000
xa = rnorm(n)
Tr = -0.5
x = xa[xa < -0.5]
tnorm.mle(x, Tr, n)
Tr = 0.5
x = xa[xa > 0.5]
tnorm.mle(x, Tr, n)
