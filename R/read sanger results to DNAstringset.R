library(sangerseqR)
library(Biostrings)
library(stringr)
library(msa)
sanger.resul.tpath <- '../'
readABI2Xstringset <- function(path){
  sanger.resul.tpath <- path
  file <- dir(sanger.resul.tpath) %>% 
    str_subset(pattern = '\\.*.ab1$') %>% 
    str_sort(numeric = T)
  seqset <- DNAStringSet(use.names=TRUE)
  for (i in seq_len(length(file))) {
    filepath <- str_c(sanger.resul.tpath,file[i])
    print(file[i])
    temp <- readsangerseq(filepath)
    seq <- as.character(temp@primarySeq)
    seqset[i] <- seq
    names(seqset)[i] <- file[i]
  }
  return(seqset)
}
seqset <- readABI2Xstringset(sanger.resul.tpath)
dsODN <- DNAString('GTTTAATTGAGTTGTCATATGTTAATAACGGTAT')
dsODN <- reverseComplement(dsODN)
dsODN <- DNAStringSet(dsODN)
names(dsODN) <- 'dsODN'
seqset <- append(seqset, dsODN, after = length(seqset))
seqset <- getAlignSeq(seqset, geneid = 'dsODN')
for (i in seq_len(12)) {
  aln <- msa(selectrow_byname(seqset, 'dsODN', names(seqset)[i]),
             method = 'ClustalOmega')
  print('*************************************************')
  print(names(seqset)[i])
  print(aln, show="complete")
}
