source('msaprint.R')
library(Biostrings)

plotFastaAln2Pdf <- function(fasta.path = filepath, 
                             title = title, 
                             type = type){
  path <- file.path(fasta.path)
  fa <- switch (type,
                DNA = readDNAStringSet(path),
                protein = readAAStringSet(path),
                RNA = readRNAStringSet(path))
  aln <- msa(fa)
  msaprintPDF(aln, title)
}
plotFastaAln2Pdf(fasta.path = '../cpf1proteins.fasta', 
                 title = 'cpf1proteins.pdf', type = 'protein')

