`plotCN` <-
function(pos, LRR, BAF, chr2plot=NULL, sampleIDs=NULL, fileNames=NULL, 
  types="genoCN", CNA=TRUE, main="", LRR.ylim=NULL, cex=0.5, plot.lowess=TRUE)
{

  if(!is.null(fileNames)){
    if(length(types) != length(fileNames)){
      stop("fileNames and types have different lengths\n")
    }
  }
  
  if(length(sampleIDs)==1){
    sampleIDs = rep(sampleIDs, length(fileNames))
  }
  
  if( any(! types %in% c("genoCN","pennCNV")) ){
    stop("unkown types\n")
  }
  
  if(any(diff(pos)<0)){
    stop("SNPs must be ordered by their positinos\n")
  }
  wNA = which(is.na(pos))
  if(!is.null(LRR)){
    wNA = union(wNA, which(is.na(LRR)))
  }
  if(!is.null(BAF)){
    wNA = union(wNA, which(is.na(BAF)))
  }

  nplot = 0
  if(!is.null(LRR)) nplot = nplot + 1
  if(!is.null(BAF)) nplot = nplot + 1
  nplot = nplot + length(fileNames)

  if(nplot>1){
    par(mfrow = c(nplot, 1))
  }
  
  if(!is.null(BAF)){
    ## plot BAF
    par(mar=c(0,4,2,2))
    plot(pos, BAF, xaxt="n", bty="n", cex.lab=1.2, cex.axis=1.1,
      main=main, cex=cex)
    abline(h = seq(0.0, 1.0, by=0.25), lty=2)
  }
  
  if(!is.null(LRR)){
    ## plot LRR
    par(mar=c(2,4,1,2))
    plot(pos, LRR, bty="n", cex.lab=1.2, cex.axis=1.1, cex=cex, ylim=LRR.ylim)
    if(plot.lowess){
      nna = which(!is.na(LRR))
      lines(lowess(LRR[nna]~pos[nna], f=1/50), lwd=2, col="red")
    }
    abline(h = 0, lty=2)
  }

  cols = c("lightblue", "darkblue", "black", "darkgreen",
   "orange", "darkred", "red", "purple", "brown")

  col.pCNV = function(states, cols){
    colp = character(length(states))
    colp[which(states==0)] = cols[3]
    colp[which(states==1)] = cols[4]
    colp[which(states==2)] = cols[2]
    colp[which(states==3)] = cols[5]
    colp[which(states==4)] = cols[7]
    colp
  }

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

  if(length(fileNames) > 0){
    for(k in 1:length(fileNames)){
      fileName = fileNames[k]
      type     = types[k]
      sampleID = sampleIDs[k]
      
      extSize  = file.info(fileName)$size

      par(mar=c(0,4,1,2))
      plot(range(pos[!is.na(pos)]), c(-0.5,4.5) , xaxt="n", yaxt="n", bty="n",
        type="n", cex.lab=1.2, cex.axis=1.2, xlab="", ylab="")
      # text(min(pos)-0.02*(max(pos) - min(pos)), 0:4, 0:4, cex=1.0)
      mtext(0:4,  side=2, line=0.5, at=0:4, cex=0.75)
      mtext(type, side=2, line=3.0, cex=0.8)

      if(extSize==0){ next }
    
      if(type=="genoCN"){
        ext = read.table(fileName, stringsAsFactors=FALSE, header=TRUE)
      }
      
      if(type=="pennCNV"){
        ext = read.table(fileName, stringsAsFactors=FALSE)
        names(ext) = c("chr", "start", "end", "state", "sample", "snp1", 
         "snp2", "score")
      }
      
      if(!is.null(sampleID)){
        ext = ext[ext$sample==sampleID,]
      }
      
      if(nrow(ext) == 0){
        start = pos[1]
        end   = pos[length(pos)]
        cn1 = 2
        col1 = cols[1]
        rect(start, cn1-0.3, end, cn1+0.3, density = NULL, angle = 45,
           col = col1, border = col1, lwd = 0.5)
        next
      }
      
      chrs = unique(ext$chr)    
      if(length(chrs)>1 & is.null(chr2plot)){
        stop("there are multiple chromosomes in input file, need to specify which
         chromosome to plot\n")
      }
      
      if(!is.null(chr2plot)){
        ext = ext[ext$chr==chr2plot,]
      }
      
      if(nrow(ext) == 0){
        start = pos[1]
        end   = pos[length(pos)]
        cn1 = 2
        col1 = cols[1]
        rect(start, cn1-0.3, end, cn1+0.3, density = NULL, angle = 45,
           col = col1, border = col1, lwd = 0.5)
        next
      }
      ext = ext[order(ext$start),]
      ww = which(pos < ext$start[1])
      if(length(ww)>0){
        start = pos[ww[1]]
        end   = pos[ww[length(ww)]]
        cn1 = 2
        col1 = cols[1]
        rect(start, cn1-0.3, end, cn1+0.3, density = NULL, angle = 45,
           col = col1, border = col1, lwd = 0.5)
      }
      
      for(i in 1:(nrow(ext)-1)){
        start = ext$start[i]
        end   = ext$end[i]
        stat  = ext$state[i]
        
        if(type=="genoCN"){
          cn1  = cn.genoCN(stat, CNA)
          col1 = cols[stat]
        }
        
        if(type=="pennCNV"){
          cn1  = stat
          col1 = col.pCNV(stat, cols)
        }
        
        rect(start, cn1-0.3, end, cn1+0.3, density = NULL, angle = 45,
             col = col1, border = col1, lwd = 0.5)
        
        ww = which(pos > end & pos < ext$start[i+1])
        if(length(ww)>0){
          start = pos[ww[1]]
          end   = pos[ww[length(ww)]]
          cn1   = 2
          col1  = cols[1]
          rect(start, cn1-0.3, end, cn1+0.3, density = NULL, angle = 45,
             col = col1, border = col1, lwd = 0.5)
        }
      }
      i = nrow(ext)
      start = ext$start[i]
      end   = ext$end[i]
      stat  = ext$state[i]
      
      if(type=="genoCN"){
        cn1  = cn.genoCN(stat, CNA)
        col1 = cols[stat]
      }
      
      if(type=="pennCNV"){
        cn1  = stat
        col1 = col.pCNV(stat, cols)
      }
      
      rect(start, cn1-0.3, end, cn1+0.3, density = NULL, angle = 45,
           col = col1, border = col1, lwd = 0.5)
      
      ww = which(pos > end)
      if(length(ww)>0){
        start = pos[ww[1]]
        end   = pos[ww[length(ww)]]
        cn1   = 2
        col1  = cols[1]
        rect(start, cn1-0.3, end, cn1+0.3, density = NULL, angle = 45,
           col = col1, border = col1, lwd = 0.5)
      }
    }
  }
}

