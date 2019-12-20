# change the fas primer file to csv ---------------------------------------

library(Biostrings)
library(stringr)
library(readr)
primerfas2csv <- function(path = path){
  primer <- readDNAStringSet(path)
  names(primer) <- str_extract(names(primer), '.*(?<=F|R)')
  primer <- as.data.frame(primer)
  path2 <- str_c(str_extract(as.character(path),'.*\\.'),'csv')
  write.csv(primer,file = path2)
}
primerfas2csv('Primers_2019_12_18.fas')

