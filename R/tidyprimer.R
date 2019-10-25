library(Biostrings)
library(stringr)
geneid <- '8085'
primerlist <- dir('../../primer/') %>% str_sort(numeric = T)
tempfile <- str_subset(primerlist, pattern = geneid)
primer <- readDNAStringSet(file.path('../../primer/', tempfile))
primer[2] <- reverseComplement(primer[2])
names(primer) <- paste(geneid,sep = '.' ,c('Fw','Rv'))
primer


