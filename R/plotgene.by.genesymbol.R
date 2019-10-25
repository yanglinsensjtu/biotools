library(Gviz)
library(biomaRt)
bm <- useMart(host="grch37.ensembl.org", 
              biomart="ENSEMBL_MART_ENSEMBL", 
              dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", 
                                    name="ENSEMBL", 
                                    symbol="TRAC", 
                                    biomart=bm)
plotTracks(biomTrack, transcriptAnnotation="symbol")

chr <- 'chr14'
f <- 22090057
t <- 23021075
bm <- useMart(host="grch37.ensembl.org", 
              biomart="ENSEMBL_MART_ENSEMBL", 
              dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", 
                                    chromosome=chr, 
                                    start=f, 
                                    end=t,
                                    name="ENSEMBL", 
                                    biomart=bm)

plotTracks(biomTrack, transcriptAnnotation="symbol")

