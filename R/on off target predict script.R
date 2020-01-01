
# Version:      0.4.0

# package install ---------------------------------------------------------

Rpackages <- c("readr",
               "tibble",
               "readxl",
               "magrittr",
               "stringr",
               "tidyr",
               "dplyr",
               "BiocManager",
               "progress")
for (i in seq_len(length(Rpackages))) {
  print(Rpackages[i])
  if (!requireNamespace(Rpackages[i], quietly = TRUE))
    install.packages(Rpackages[i])
}




bioconductorpackage <- c("Biostrings",
                         "GenomicRanges",
                         "TxDb.Hsapiens.UCSC.hg38.knownGene",
                         "VariantAnnotation",
                         "BSgenome.Hsapiens.UCSC.hg38",
                         "Gviz",
                         "biomaRt",
                         "AnnotationHub")
for (i in seq_len(length(bioconductorpackage))) {
  print(bioconductorpackage[i])
  if (!requireNamespace(bioconductorpackage[i], quietly = TRUE))
    BiocManager::install(bioconductorpackage[i])
}

# read the raw data into the R --------------------------------------------
library(readr)
library(tibble)
library(readxl)
library(magrittr)

off.spotter<- read.delim('offspotterhg38.txt',
                         header = F,
                         stringsAsFactors = F) %>%
  as_tibble()

cas.offfinder <- read.delim('casoffinderhg38.txt',
                            header = T,
                            stringsAsFactors = F) %>%
  as_tibble()

cosmid <- read_xlsx('cosmidhg38.xlsx')

# transform the data into granges obj -------------------------------------

library(GenomicRanges)
library(stringr)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqinfo <- Seqinfo(genome = 'hg38')
off.spotter$V1 <- str_c('chr',off.spotter$V1)
names(off.spotter)[5:6] <- c('query','hit')
off.spotter.gr <- GRanges(seqnames = off.spotter$V1,
                          strand = off.spotter$V2,
                          ranges = IRanges(start = off.spotter$V3,
                                           end = off.spotter$V4),
                          seqinfo = seqinfo,
                          dplyr::select(off.spotter,
                                        query,
                                        hit))


cas.offfinder.gr <- GRanges(seqnames = cas.offfinder$Chromosome,
                            strand = cas.offfinder$Direction,
                            ranges = IRanges(start = cas.offfinder$Position,
                                             width = 23),
                            seqinfo = seqinfo,
                            dplyr::select(cas.offfinder,
                                          X.Bulge.type,
                                          crRNA,
                                          DNA,
                                          Mismatches,
                                          Bulge.Size))

cosmid <- separate(cosmid, `Chr Position`, sep = ':', into = c('chromosome', 'position'))
cosmid$chromosome <- stringr::str_replace_all(cosmid$chromosome, 'Chr', 'chr')
cosmid <- separate(cosmid, position, sep = '-', into = c('position.start', 'position.end'))
names(cosmid) <- str_replace_all(names(cosmid), ' ', '.')
cosmid.gr <- GRanges(seqnames = cosmid$chromosome,
                     strand = str_sub(cosmid$Strand,1,1),
                     ranges = IRanges(start = as.integer(cosmid$position.start),
                                      end = as.integer(cosmid$position.end)),
                     seqinfo = seqinfo,
                     dplyr::select(cosmid,
                                   Score,
                                   Query.tag,
                                   Search.result,
                                   Query.type,
                                   Mismatch,
                                   Ends.with.RG,
                                   Cut.site))


cc <- unique(c(cosmid.gr, off.spotter.gr,cas.offfinder.gr))

# Annotation --------------------------------------------------------------

library(AnnotationHub)
library(GenomicFeatures)
library(GenomicRanges)
ah = AnnotationHub()
GRCh38.96 <- ah[['AH69461']]
#GRCh38.96.chr <- ah[['AH69459']]
org.NCBI <- ah[['AH70572']]
seqlevels(GRCh38.96) <- paste0('chr', seqlevels(GRCh38.96))
genome(GRCh38.96) <- 'hg38'



return.overlap <- function(c1, c2, minoverlap=0L, ignore.strand = FALSE)
{
  if (unique(genome(c1)) == unique(genome(c2)))
  {
    if (as.character(c1@seqnames) == as.character(c2@seqnames))
    {
      if (as.character(c1@strand) == as.character(c2@strand))
      {
        tmpOver <- subsetByOverlaps(c1@ranges, c2@ranges, minoverlap = minoverlap)
        j <- isEmpty(tmpOver)
        if(!j)
        {
          return(c1)
        }else
        {
          return(NA)
        }
      }else if (ignore.strand)
      {
        tmpOver.ignore <- subsetByOverlaps(c1@ranges, c2@ranges, minoverlap = minoverlap)
        j.ignore <- isEmpty(tmpOver.ignore)
        if (!j.ignore)
        {
          return(c1)
        }else
        {
          return(NA)
        }
      }else {
        #print("There is no overlap")
        return(NA)
      }
    }else
    {
      #print("The chromosome was not the same")
      return(NA)
    }
  }else
  {
    #print("The genome was not the same")
    return(NA)
  }
}

library(progress)
mapsgRNA2CDS <- function(predict.granges = granges)
{
  var <- predict.granges
  gr <- GRanges()
  pb <- progress_bar$new(
    format = "完成百分比 [:bar] :percent 执行时间 :elapsed",
    total = length(var), clear = FALSE, width= 60)
  for (i in seq_len(length(var))) {
    pb$tick()
    tmp <- var[i]
    map <- subsetByOverlaps(GRCh38.96, tmp, ignore.strand=T)
    if("CDS" %in% map$type)
    {
      tmp$type = "CDS"
      tmp$gene_name = unique(map[map$type =="CDS"]$gene_name)
      tmp$gene_id = unique(map[map$type =="CDS"]$gene_id)
      #此处添加判断，如果已在集合中，丢弃，如果不在集合中，则加入集合中。
      j <- FALSE
      for (i in seq_along(gr))
      {
        if(!is.na(return.overlap(gr[i], tmp, ignore.strand = T, minoverlap = 17)))
        {
          j <- TRUE
          break()
        }
      }
      if(!j)
      {
        gr <- append(gr, tmp, after = length(gr))
      }

    }
  }
  gru <- unique(gr)
  return(gru)
}





cc.ann <- mapsgRNA2CDS(cc)
keys <- cc.ann$gene_name
GENESy2ID <- select(org.NCBI,
                    keys = keys,
                    columns = c("ENTREZID","SYMBOL"),
                    keytype = "SYMBOL")
cc.ann$GENEID <- GENESy2ID$ENTREZID
cco.ann <- cc.ann

# read the sgRNA sequence into R ------------------------------------------
library(Biostrings)
sgRNA <- readDNAStringSet('sgRNAsequence.txt')

# agRNA alignment function ------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg38)
BS.hg38 <- BSgenome.Hsapiens.UCSC.hg38
# sgRNA aliganment --------------------------------------------------------
sgRNA.alignment <- function(grange = Grange.obj,
                            sgRNA = targetgene){
  off.target.seq <- getSeq(BS.hg38, unique(grange))
  if (is.null(grange$GENEID)) {
    print('There is no geneid')
    for (i in seq_len(length(unique(grange)))) {
      print(i)
      print(pairwiseAlignment(sgRNA,
                              off.target.seq[i],
                              gapOpening = 0,
                              gapExtension = 1))
    }
  }else{
    for (i in seq_len(length(unique(grange)))) {
      print(i)

      if (is.na(grange$GENEID[i])) {
        print('The geneid is NA')
      }else{
        print(grange$GENEID[i])
      }

      print(pairwiseAlignment(sgRNA,
                              off.target.seq[i],
                              gapOpening = 0,
                              gapExtension = 1))
    }
  }
}

sgRNA.alignment(cco.ann,sgRNA)

library(Biostrings)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg38)
BS.hg38 <- BSgenome.Hsapiens.UCSC.hg38

# write the alignment to txt fucntion -------------------------------------


write.alignment2txt <- function(granges.ann = granges.ann,
                                sgRNA = sgRNA,
                                path = path){
  print(granges.ann$GENEID)
  print(granges.ann$GENEID)
  off.target.seq <- getSeq(BS.hg38, granges.ann)
  write_lines('****The sgRNA alignment with the off target sites***',
              path = path)
  for (i in seq_len(length(granges.ann))) {
    cat(paste(i, ' '))
    paln <- pairwiseAlignment(sgRNA,
                              off.target.seq[i],
                              gapOpening = 0,
                              gapExtension = 1)

    write_lines('------------------------',
                path = path,
                append = T)
    write_lines(paste('--',i,'--'),
                path = path,
                append = T)
    write_lines(names(sgRNA),
                path = path,
                append = T)
    write_lines(paste(granges.ann$GENEID[i],
                      granges.ann$gene_name[i],
                      as.character(granges(granges.ann)[i])),
                path = path,
                append = T)
    write_lines('.....',
                path = path,
                append = T)
    write_lines(alignedPattern(paln),
                path = path,
                append = T)
    write_lines(alignedSubject(paln),
                path = path,
                append = T)

  }
}

write.alignment2txt(cco.ann,
                    sgRNA,
                    path = 'align sgRNA target.txt')

# write the 600bp offtarget site into txt ---------------------------------
library(GenomicRanges)
# change granges obj ranges -----------------------------------------------

changeIRanges <- function(granges.obj = granges.obj,
                          upstream = integer,
                          width = integer, ...)
{
  if (as.character(granges.obj@strand) == "-")
  {
    GRanges(seqnames = granges.obj@seqnames,
            IRanges(end = (end(granges.obj@ranges) + upstream),
                    start = (end(granges.obj@ranges) + upstream - width + 1),
                    width = width),
            strand = granges.obj@strand,
            seqinfo = granges.obj@seqinfo,
            mcols(granges.obj))
  } else
  {
    GRanges(seqnames = granges.obj@seqnames,
            IRanges(start = start(granges.obj@ranges) - upstream,
                    end = start(granges.obj@ranges) - upstream + width -1,
                    width = width),
            strand = granges.obj@strand,
            seqinfo = granges.obj@seqinfo,
            mcols(granges.obj))
  }
}

cco.ann.extend <- changeIRanges(granges.obj = cco.ann,
                                upstream = 300,
                                width = 600)

cco.ann.extend.seq <- getSeq(BS.hg38, cco.ann.extend)

names(cco.ann.extend.seq) <- paste0(c(1:length(cco.ann.extend$GENEID)),
                                    "-",
                                    cco.ann.extend$GENEID,
                                    "-",
                                    cco.ann.extend$gene_name)

writeXStringSet(cco.ann.extend.seq,
                filepath = 'predict off target genes sequences on the genome.txt')


# map the offtarget to the genome -----------------------------------------

library(Gviz)
library(biomaRt)
library(stringr)

plotGviz <- function(coding.obj, i, bounds, plottitle = '', folder){
  tempfolder <- file.path(folder)

  if (!file.exists(tempfolder)) {
    print(paste('The directory does not exist', folder, 'file folder'))
    print('Creating folder!')
    dir.create(tempfolder)
  }
  chr <- as.character(coding.obj@seqnames[i])
  strand <- as.character(coding.obj@strand[i])
  start <- IRanges::start(coding.obj)[i]
  end <- IRanges::end(coding.obj)[i]
  filetitle <- str_c(i,chr,'-',start,'-',end,'.jpg')
  print(filetitle)
  # annotation track creation -----------------------------------------------

  anntrack  <- AnnotationTrack(start = start,
                               end = end,
                               strand = strand,
                               chromosome = chr,
                               genome ="hg38",
                               name ="CRISPR")

  # biomartGeneRegiontrack creation -----------------------------------------

  f <- as.integer(start - bounds)
  t <- as.integer(start + bounds)
  bm <- useMart(host="www.ensembl.org",
                biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")
  biomTrack <- BiomartGeneRegionTrack(genome="hg38",
                                      chromosome=chr,
                                      start=f,
                                      end=t,
                                      name="ENSEMBL",
                                      biomart=bm)

  # Ideogramtrack creation --------------------------------------------------

  itrack <- IdeogramTrack(genome = 'hg38',chromosome = chr)

  # genome axis track creation ----------------------------------------------

  gatrack <- GenomeAxisTrack(distFromAxis = 15,labelPos="below")

  # save the image ----------------------------------------------------------
  oldpath <- getwd()
  setwd(dir = tempfolder)
  sz <- c(1,1,4,1)
  ls <- list(itrack,gatrack,biomTrack,anntrack)
  mn <- str_c(i,'-',chr,'...',plottitle)
  jpeg(filetitle, width = 5000, height =3090,res = 720)
  plotTracks(trackList = ls,
             transcriptAnnotation="symbol",
             sizes = sz,
             from = f,
             to = t,
             main = mn,
             cex.main = 0.8)
  dev.off()

  setwd(dir = oldpath)
}

coding.obj <- cco.ann
GvizVar <- length(coding.obj)
pbGivz <- progress_bar$new(
  format = "完成百分比 [:bar] :percent 执行时间 :elapsed",
  total = GvizVar, clear = FALSE, width= 60)
for (i in seq_along(coding.obj)) {
  pbGivz$tick()
  chr <- as.character(coding.obj@seqnames[i])
  strand <- as.character(coding.obj@strand[i])
  start <- IRanges::start(coding.obj)[i]
  end <- IRanges::end(coding.obj)[i]
  filetitle <- str_c(i,chr,'-',start,'-',end,'.jpg')
  if (file.exists(paste0('on off target site map to genome/',filetitle))) {
    print(paste(filetitle,'exists'))
    next()
  }else{
    try(plotGviz(coding.obj = coding.obj,
             i=i,
             bounds = 10000,
             plottitle = str_c(coding.obj$gene_name[i],
                               '-',
                               coding.obj$GENEID[i],
                               '-',
                               coding.obj$gene_id[i]),
             folder = 'on off target site map to genome'))
  }

}

