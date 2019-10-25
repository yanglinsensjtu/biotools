library(sangerseqR)
library(Biostrings)
library(stringr)

source('sangerseqquility.R')
sanger.resul.tpath <- '../sanger seq results/'

filelist <- dir(sanger.resul.tpath) %>% 
  str_subset(pattern = '\\.*.ab1$') %>% stringr::str_sort(numeric = T)

nonspecificfolder <- file.path(sanger.resul.tpath, 'nonspecific')

if (file.exists(nonspecificfolder)) {
  print('The nonspecific folder already existed')
}else{
  dir.create(nonspecificfolder)
}


for (i in seq_len(length(filelist))) {
  filepath <- str_c(sanger.resul.tpath,filelist[i])
  temp <- readsangerseq(filepath)
  temp <- score.sanger(temp)
  score <- mean(temp$score)
  if (score > 2) {
    next()
  }else{
    file.copy(filepath, nonspecificfolder)
    file.remove(filepath)
  }
  
}
