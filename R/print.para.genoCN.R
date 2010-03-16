`print.para.genoCN` <-
function(x, digits=3, ...)
{
  if(x$CNA){
    message("\nParameters for CNA studies\n\n")   
  }else{
    message("\nParameters for CNV studies\n\n")   
  }
  
  for(i in 1:length(x)){
  
    if(is.character(x[[i]])){
      message("------------------------------------------\n")
      message(names(x)[i], ": ", sep="")
      message(x[[i]], "\n")
    }else if(is.vector(x[[i]])){
      message("------------------------------------------\n")
      message(names(x)[i], ": ", sep="")
      message(round(x[[i]], digits), "\n")
    }else{
      message("------------------------------------------\n")
      message(names(x)[i], ":\n", sep="")
      print(round(x[[i]], digits))
    }
  }
  
  
  message("------------------------------------------\n")
}
