library(Biostrings)
library(msa)

getProteinSeq <- function(uniprotProteinID){
  url <- paste0('https://www.uniprot.org/uniprot/',uniprotProteinID, '.fasta')
  print(url)
  fastaName <- paste0(uniprotProteinID,'.fasta')
  download.file(url, fastaName)
  proteinSeq <- readDNAStringSet(fastaName)
  filepath <- paste('The file', fastaName, 'is save in',getwd(), '!!!')
  return(proteinSeq)
}

source('msaprint.R')
B2M <- readDNAStringSet('../B2M.fasta')
P61769 <- getProteinSeq('P61769')
u <- append(P61769 ,B2M)
u
aln <- msa(u,
           method = 'ClustalW'
             )
print(aln, show = 'complete')
msaprintPDF(align.obj = aln, title = 'trbc2.pdf')
