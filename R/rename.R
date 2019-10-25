library(stringr)
sanger.resul.tpath <- '../sanger seq results/1270/'
filelist <- dir(sanger.resul.tpath) %>% 
  str_subset(pattern = '\\.*.ab1$')
geneid <- '1270'
pat <- str_c(geneid,'\\.*')
file <- str_subset(filelist, pattern = pat)
file <- str_sort(file, numeric = T)
file
for (i in seq_len(length(file))) {
  filepath <- str_c(sanger.resul.tpath,file[i])
  rename <- str_replace(file[i],'1270','1043')
  filepath.rename <- str_c(sanger.resul.tpath, rename)
  file.rename(from = filepath,to = filepath.rename)
}

