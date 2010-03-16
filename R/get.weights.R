# ----------------------------------------------------------------------
# return a matrix of size L*M, where L is the number of SNPs and
# M is the number of genotype classes. M=36 for CNA studies and 
# M=21 for CNV studies
# ----------------------------------------------------------------------

get.weights <- function(p.Bs, CNA, contamination, normalGtp, geno.error)
{
  genoP = 1-geno.error
  
  p.As = 1 - p.Bs
  L    = length(p.Bs)
  if(CNA){ M = 36 } else{ M = 21 }
  
  ws   = matrix(0, nrow=L, ncol=M)

  ##
  #  the columns of matrix ws
  #  HMM state  Genotype state
  #  1          1 -> AA      2 -> AB        3 -> BB
  #  2          4 -> AA      5 -> (AA,AB)   6 -> (BB,AB)   7 -> BB
  #  3          8 -> NULL
  #  4          9 -> A      10 -> (A,AB)   11 -> (B,AB)   12 -> B
  #  5          13 -> AAA   14 -> AAB      15 -> ABB      16 -> BBB
  #  6/CNV      17 -> AAAA  18 -> AAAB     19 -> AABB     20 -> ABBB 21 -> BBBB
  #  6/CNA      22 -> AAA   23 -> (AAA,AB) 24 -> (BBB,AB) 25 -> BBB
  #  7          26 -> AAAA  27 -> AABB     28 -> BBBB
  #  8          29 -> AAAA  30 -> (AAAA,AB)31 -> (BBBB,AB)32 -> BBBB
  #  9          33 -> AAAA  34 -> AAAB     35 -> ABBB     36 -> BBBB
  ##

  # state 3 (NULL)
  ws[,8] = rep(1,L)
  
  if(!CNA){
    ## ----------------------------------------------
    ## IF it is for CNV studies
    ## ----------------------------------------------

    # state 1 (AA, AB, BB)
    ws[,1] = p.As*p.As
    ws[,2] = 2*p.As*p.Bs
    ws[,3] = p.Bs*p.Bs
    
    # state 2 and 4
    ws[,c(4,9)]  = p.As
    ws[,c(7,12)] = p.Bs

    # state 5 (AAA, AAB, ABB, BBB)
    ws[,13] = dbinom(0,3,p.Bs)
    ws[,14] = dbinom(1,3,p.Bs)
    ws[,15] = dbinom(2,3,p.Bs)
    ws[,16] = dbinom(3,3,p.Bs)

    # state 6 (AAAA, AAAB, AABB, ABBB, BBBB)
    ws[,17] = dbinom(0,4,p.Bs)
    ws[,18] = dbinom(1,4,p.Bs)
    ws[,19] = dbinom(2,4,p.Bs)
    ws[,20] = dbinom(3,4,p.Bs)
    ws[,21] = dbinom(4,4,p.Bs)
    
  }else{
    ## ----------------------------------------------
    ## IF it is for CNA sutides
    ## ----------------------------------------------

    Q = which(normalGtp == -1)
  
    if(length(Q) > 0){
      p.Bs.Q = p.Bs[Q]
      p.As.Q = 1 - p.Bs.Q
            
      # state 1 (AA, AB, BB)
      ws[Q,1] = p.As.Q*p.As.Q
      ws[Q,2] = 2*p.As.Q*p.Bs.Q
      ws[Q,3] = p.Bs.Q*p.Bs.Q
  
      # state 2 (AA, (AA,AB), (BB,AB), BB)
      # state 4 (A, (A,AB), (B,AB), B)
      if(contamination){
        ws[Q,c(4,9)]       = ws[Q,1]
        ws[Q,c(5:6,10:11)] = ws[Q,2]/2
        ws[Q,c(7,12)]      = ws[Q,3]
      }else{
        ws[Q,c(4,9)]  = p.As.Q
        ws[Q,c(7,12)] = p.Bs.Q
      }
      
      # state 5 (AAA, AAB, ABB, BBB)
      ws[Q,13]    = ws[Q,1]
      ws[Q,14:15] = ws[Q,2]/2
      ws[Q,16]    = ws[Q,3]
  
      # state 6 (AAA, (AAA, AB), (AB, BBB), BBB)
      if(contamination){
        ws[Q,22]    = ws[Q,1]
        ws[Q,23:24] = ws[Q,2]/2
        ws[Q,25]    = ws[Q,3]
      }else{
        ws[Q,22] = p.As.Q
        ws[Q,25] = p.Bs.Q
      }
      
      # state 7 (AAAA, AABB, BBBB)
      ws[Q,26]    = ws[Q,1]
      ws[Q,27]    = ws[Q,2]
      ws[Q,28]    = ws[Q,3]
      
      # state 8 (AAAA, (AAAA,AB), (BBBB,AB), BBBB)
      if(contamination){
        ws[Q,29]    = ws[Q,1]
        ws[Q,30:31] = ws[Q,2]/2
        ws[Q,32]    = ws[Q,3]
      }else{
        ws[Q,29]    = p.As.Q
        ws[Q,32]    = p.Bs.Q      
      }
      
      # state 9 (AAAA, AAAB, ABBB, BBBB)
      ws[Q,33]    = ws[Q,1]
      ws[Q,34:35] = ws[Q,2]/2
      ws[Q,36]    = ws[Q,3]
    }
    
    ## ----------------------------------------------
    # normal genotype is "AA"
    ## ----------------------------------------------
  
    Q = which(normalGtp==0)
    if(length(Q) > 0){
      ## weight for genotype in state 1,2,4-9
      ## given normal genotype data
      ws[Q,c(1,4,9,13,22,26,29,33)] = genoP
      
      ## state 1
      ws[Q,2:3] = geno.error/2
  
      ## state 2/4
      if(contamination){
        ws[Q,c(5:7,10:12)] = geno.error/3
      }else{
        ws[Q,c(7,12)] = geno.error
      }
      
      ## state 5
      ws[Q,14:16] = geno.error/3
      
      ## state 6
      if(contamination){
        ws[Q,23:25] = geno.error/3
      }else{
        ws[Q,25] = geno.error
      }
  
      ## state 7
      ws[Q,27:28] = geno.error/2
      
      ## state 8
      if(contamination){
        ws[Q,30:32] = geno.error/3
      }else{
        ws[Q,32] = geno.error
      }
      
      ## state 9
      ws[Q,34:36] = geno.error/3
    }
     
    # normal genotype is "AB"
    Q = which(normalGtp==1)
    if(length(Q) > 0){
      ## state 1
      ws[Q,2] = genoP
      ws[Q,c(1,3)] = geno.error/2
  
      ## state 2/4
      if(contamination){
        ws[Q,c(5,6,10,11)] = genoP/2
        ws[Q,c(4,7,9,12)]  = geno.error/2
      }else{
        ws[Q,c(4,7,9,12)] = 0.5
      }
      
      ## state 5
      ws[Q,14:15] = genoP/2
      ws[Q,c(13,16)] = geno.error/2
      
      ## state 6
      if(contamination){
        ws[Q,23:24] = genoP/2
        ws[Q,c(22,25)] = geno.error/2
      }else{
        ws[Q,c(22,25)] = 0.5      
      }
      
      ## state 7
      ws[Q,27] = genoP
      ws[Q,c(26,28)] = geno.error/2
      
      ## state 8
      if(contamination){
        ws[Q,30:31] = genoP/2
        ws[Q,c(29,32)] = geno.error/2
      }else{
        ws[Q,c(29,32)] = 0.5      
      }
      
      ## state 9
      ws[Q,34:35] = genoP/2
      ws[Q,c(33,36)] = geno.error/2
    }
  
    # normal genotype is "BB"
    Q = which(normalGtp==2)
    if(length(Q) > 0){
      ws[Q,c(3,7,12,16,25,28,32,36)] = genoP
      
      ## state 1
      ws[Q,1:2] = geno.error/2
  
      ## state 2/4
      if(contamination){
        ws[Q,c(4:6,9:11)] = geno.error/3
      }else{
        ws[Q,c(4,9)] = geno.error
      }
      
      ## state 5
      ws[Q,13:15] = geno.error/3
      
      ## state 6
      if(contamination){
        ws[Q,22:24] = geno.error/3
      }else{
        ws[Q,22] = geno.error
      }
      
      ## state 7
      ws[Q,26:27] = geno.error/2

      ## state 8
      if(contamination){
        ws[Q,29:31] = geno.error/3
      }else{
        ws[Q,29] = geno.error
      }
      
      ## state 9
      ws[Q,33:35] = geno.error/3
    }
  }
  
  ws
}
