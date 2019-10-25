library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
BS19 <- BSgenome.Hsapiens.UCSC.hg19
columns(txdb)
keytypes(txdb)
geneid <- '567'
chrv <- c( "CDSCHROM" ,
           "CDSEND"  ,
           "CDSID",
           "CDSNAME",
           "CDSSTART",
           "CDSSTRAND", 
           "EXONCHROM",
           "EXONEND",
           "EXONID",
           "EXONNAME",
           "EXONRANK",
           "EXONSTART", 
           "EXONSTRAND",
           "GENEID" )
B2M <- select(txdb, keys = geneid, columns = chrv, keytype = 'GENEID')
exons(B2M)
B2Mexons <- exons(txdb, filter = list(gene_id= '567'))
B2Mexonsseq <- getSeq(BS19,B2Mexons)
names(B2Mexonsseq) <- B2Mexons$exon_id
writeXStringSet(B2Mexonsseq, 'B2Mexon.fasta')
