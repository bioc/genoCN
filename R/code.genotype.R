`code.genotype` <-
function(v){
  if(!is.character(v)){
    stop("v must be a characer vector\n")
  }
  
  gs = unique(v)
  gs = gs[!is.na(gs)]

  alleles = unique(unlist(strsplit(gs, split="")))

  if(length(alleles) == 2){
    a1 = alleles[1]
    a2 = alleles[2]
  }else if(length(alleles) == 1){
    a1 = alleles[1]
    ## warning("only one unique allele found \n")
    return(rep(0, length(v)))
  }else{
    stop("the number of unique alleles is neither 1 nor 2\n")
  }

  g1 = paste(a1, a1, sep="")
  g3 = paste(a2, a2, sep="")
  g21 = paste(a1, a2, sep="")
  g22 = paste(a2, a1, sep="")

  which.g1 = which(v==g1)
  which.g2 = which(v==g21 | v==g22)
  which.g3 = which(v==g3)
  which.na = which(is.na(v))

  len.g1 = length(which.g1)
  len.g2 = length(which.g2)
  len.g3 = length(which.g3)
  len.na = length(which.na)

  if(len.g1 + len.g2 + len.g3 + len.na != length(v)){
    stop("wait, I miss something here..:[\n")
  }

  ndrow = rep(1, length(v))
  ndrow[which.na] = NA
  if(len.g1 >= len.g3){
    ndrow[which.g1] = 0
    ndrow[which.g3] = 2
  }else{
    ndrow[which.g1] = 2
    ndrow[which.g3] = 0
  }

  ndrow
}

