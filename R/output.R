
# ----------------------------------------------------------------------
#  write out the results
# ----------------------------------------------------------------------

output <- function(object, snpNames, chr1, pos, sampleID, header, 
  output.loh, CNA, tag="genoCN_output", outputSeg=TRUE, outputSNP=1, 
  seg.nSNP=3)
{
  
  if(class(object) != "postP.genoCN"){
    stop("the parameter must be of class 'postP.genoCN'\n")
  }
   
  # if there is no CNV probes
  if(header){
    if(outputSNP==1 || outputSNP==2){
      colNames = c("name", "state", "stateP", "CN", "CNP", "genotype", "genoP")
      cat(colNames, file = sprintf("%s_SNP.txt", tag), sep = "\t")
      cat("\n", file = sprintf("%s_SNP.txt", tag), append=TRUE)
    }
    
    if(outputSNP==3){
      colNames = c("name", "cn0", "cn1", "cn2", "cn3", "cn4", "AA", "AB", "BB")
      colNames = c(colNames, c( "Null", "A", "B", "AAA", "AAB", "ABB"))
      colNames = c(colNames, c("BBB", "AAAA", "AAAB", "AABB", "ABBB", "BBBB"))
      cat(colNames, file = sprintf("%s_SNP.txt", tag), sep = "\t")
      cat("\n", file = sprintf("%s_SNP.txt", tag), append=TRUE)
    }
    
    if(outputSeg){
      colNames = c("chr", "start", "end", "state", "cn", "sample", "snp1", "snp2")
      colNames = c(colNames, c("score", "n"))
      cat(colNames, file = sprintf("%s_segment.txt", tag), sep = "\t")
      cat("\n", file = sprintf("%s_segment.txt", tag), append=TRUE)
    }
  }
  
  HMMstate = object$state.max$state
  
  if(output.loh){
    wkp = which(HMMstate != 1)
  }else{
    wkp = which(HMMstate != 1 & HMMstate != 2)
  }

  if(length(wkp)==0){ return(0) }
  
  if(outputSNP){
    if(outputSNP==1){
      out = data.frame(name=snpNames[wkp], object$state.max[wkp,], 
       object$cn.max[wkp,], object$geno.max[wkp,])
    }else if(outputSNP==2){
      out = data.frame(name=snpNames, object$state.max, 
       object$cn.max, object$geno.max)
    }else if(outputSNP==3){
      out = data.frame(name=snpNames, object$cn, object$genotype)
    }else{
      stop("outputSNP must be 0, 1, 2, or 3\n")
    }
    
    if(outputSNP==1 || outputSNP==2){
      # column 3, 5, 7 are "stateP", "CNP", and "genoP".
      for(i in c(3,5,7)){ out[[i]] = round(out[[i]], 7) }
    }else if(outputSNP==3){
      for(i in 2:ncol(out)){ out[[i]] = round(out[[i]], 7) }
    }
    
    write.table(out, file = sprintf("%s_SNP.txt", tag), quote = FALSE, 
      sep = "\t", row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  
  # ----------------------------------------------------------------------
  # wkp is indimessageor of the SNPs with copy number is not 2
  # the segments boundary (bd) are defined as those points where 
  # wkp are two un-adjacnet SNPs or if the ajacent SNPs in wkp  
  # have different states 
  # ----------------------------------------------------------------------
  if(outputSeg){
    pos.wkp = pos[wkp]
    snp.wkp = snpNames[wkp]
    sat.wkp = HMMstate[wkp]
    
    sat.bd  = which(sat.wkp[-1] != sat.wkp[-length(sat.wkp)])
    bd      = sort(union(which(diff(wkp) > 1), sat.bd))
    starts  = c(1, bd+1)
    ends    = c(bd, length(wkp))
    
    nn        = length(starts)
    seg.n     = ends - starts + 1
    seg.start = pos.wkp[starts]
    seg.end   = pos.wkp[ends]
    seg.sat   = sat.wkp[starts]
    seg.snp1  = snp.wkp[starts]
    seg.snp2  = snp.wkp[ends]
    statepP   = object$state.max$maxpP[wkp]
    
    cn.genoCN = function(states, CNA){
      cn1 = states
      cn1[which(states==1 | states==2)] = 2
      cn1[which(states==3)] = 0
      cn1[which(states==4)] = 1
      if(CNA){
        cn1[which(states %in% 5:6)] = 3
        cn1[which(states %in% 7:9)] = 4
      }else{
        cn1[which(states==5)] = 3
        cn1[which(states==6)] = 4
      }
      cn1
    }
    
    seg.cn  = cn.genoCN(seg.sat, CNA)
    
    seg.score = numeric(nn)
    for(i in 1:nn){
      seg.score[i] = sum(statepP[starts[i]:ends[i]])
    }
    seg.score = round(seg.score,5)
    
    seg = data.frame(chr=rep(chr1, length(seg.start)), 
      start=seg.start, end=seg.end, state=seg.sat, 
      cn=seg.cn, sample=rep(sampleID,nn), snp1=seg.snp1, 
      snp2=seg.snp2, score=seg.score, n=seg.n)
    
    seg = seg[seg$n >= seg.nSNP,]
    
    write.table(seg, file = sprintf("%s_segment.txt", tag), 
                quote = FALSE, sep = "\t", row.names = FALSE, 
                col.names = FALSE, append=TRUE)
  }
}

