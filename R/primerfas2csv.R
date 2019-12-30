# change the fas primer file to csv ---------------------------------------

library(Biostrings)
library(stringr)
library(readr)
primerfas2csv <- function(path = path, prefix = NA){
  prefix <- prefix
  primer <- readDNAStringSet(path)
  names(primer) <- str_extract(names(primer), '.*(?<=F|R)')
  if(!is.na(prefix))
  {
    names(primer) <- str_c(prefix,names(primer))
  }
  primer <- as.data.frame(primer)
  path2 <- str_c(str_extract(as.character(path),'.*\\.'),'csv')
  write.csv(primer,file = path2)
}
primerfas2csv('Primers_2019_12_30.fas', prefix = 'C2')
