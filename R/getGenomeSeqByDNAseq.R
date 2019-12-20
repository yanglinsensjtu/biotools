library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19



changeIRanges <- function(granges.obj = granges.obj, upstream = integer, width = integer){
  GRanges(seqnames = granges.obj@seqnames,
          IRanges(start = start(granges.obj@ranges) - upstream,
                  width = width),
          strand = granges.obj@strand,
          seqinfo = granges.obj@seqinfo,
          mcols(granges.obj))
}

getSeqFgenome <- function(seqstr = sequence,title = title, ...){
  seqstr <- seqstr
  seq <- DNAString(seqstr)
  chr.name <- seqnames(BS.hg19)
  #chr.name <- sample(chr.name, length(chr.name), replace = F)
  Flank.seq <- NA
  for (i in seq_len(length(chr.name))) {
    chr.seq <- getSeq(BS.hg19, chr.name[i])
    seqr <- reverseComplement(seq)
    mtp <- matchPattern(pattern = seq, chr.seq)
    mtpr <- matchPattern(pattern = seqr, chr.seq)
    if(length(mtp)!=0){
      ranges <- mtp@ranges
      strand <- '+'
      location <- GRanges(seqnames = chr.name[i],
                          ranges = ranges,
                          strand = strand)
      location.c <- changeIRanges(location, ...)
      Flank.seq <- getSeq(BS.hg19, location.c)
      names(Flank.seq) <- paste0(title,
                                 ' ',
                                 location.c@seqnames,
                                 ' ',
                                 location.c@strand,
                                 ':',
                                 ' ',
                                 location.c@ranges@start,
                                 '-',
                                 location.c@ranges@start)
      print(paste0('The query sequence is sense strand according to the Genome ',chr.name[i]))
      break()
    }else if(length(mtpr)!=0){
      rangesr <- mtpr@ranges
      strand <- '-'
      location <- GRanges(seqnames = chr.name[i],
                          ranges = ranges,
                          strand = strand)

      location.c <- changeIRanges(location, ...)
      Flank.seq <- getSeq(BS.hg19, location.c)
      names(Flank.seq) <- paste0(title,
                                 ' ',
                                 location.c@seqnames,
                                 ' ',
                                 location.c@strand,
                                 ':',
                                 ' ',
                                 location.c@ranges@start,
                                 '-',
                                 location.c@ranges@start)
      print(paste0('The query sequence is anti-sense strand according to the BS.hg19 Genome ',chr.name[i]))
      break()
    }else{
      print(paste0('Not found any location match the sequence in ',chr.name[i]))
    }
  }

  if (!is.na(Flank.seq)) {

    writeXStringSet(Flank.seq, filepath = paste0(title,'.fasta'))
    return(Flank.seq)
  }else{
    return(NA)
  }

}

getSeqFgenome('', title = '',upstream = 100, width = 200)





