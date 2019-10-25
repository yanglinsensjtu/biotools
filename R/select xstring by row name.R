#select row of xstring object by name
library(stringr)
library(Biostrings)

selectrow_byname <- function(x = xstring, 
                             xstringname, 
                             ...){
  temp <- x[c(xstringname, ...)]
  return(temp)
}

