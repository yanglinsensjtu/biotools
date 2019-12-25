library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19

changeIRanges <- function(granges.obj = granges.obj, upstream = integer, width = integer){
  if (as.character(granges.obj@strand) == '-') {
    GRanges(seqnames = granges.obj@seqnames,
            IRanges(end = (end(granges.obj@ranges) + upstream),
                    start = (end(granges.obj@ranges) + upstream - width + 1)),
            strand = granges.obj@strand,
            seqinfo = granges.obj@seqinfo,
            mcols(granges.obj))
  }else{
    GRanges(seqnames = granges.obj@seqnames,
            IRanges(start = start(granges.obj@ranges) - upstream ,
                    width = width),
            strand = granges.obj@strand,
            seqinfo = granges.obj@seqinfo,
            mcols(granges.obj))
  }
}
getSeqFgenome <- function(seqstr = sequence, 
                          title = NA, 
                          chrMatch = NA, 
                          upstream = 0, 
                          width = NA,
                          ...){
  
  seqstr <- seqstr
  if (is.na(width)) {
    width <- nchar(seqstr)
  }
  seq <- DNAString(seqstr)
  Flank.seq <- NA
  if (is.na(chrMatch)) {
    chr.name <- seqnames(BS.hg19)
    chr.name <- sample(chr.name, length(chr.name), replace = F)
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
        location.c <- changeIRanges(location, upstream, width)
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
                                   (location.c@ranges@width + location.c@ranges@start -1 ))
        print(paste0('The query sequence is sense strand according to the Genome ',chr.name[i]))
        #break()
      }else if(length(mtpr)!=0){
        rangesr <- mtpr@ranges
        strand <- '-'
        location <- GRanges(seqnames = chr.name[i],
                            ranges = rangesr,
                            strand = strand)
        
        location.c <- changeIRanges(location, upstream, width)
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
                                   (location.c@ranges@width + location.c@ranges@start ))
        print(paste0('The query sequence is anti-sense strand according to the BS.hg19 Genome ',chr.name[i]))
        #break()
      }else{
        print(paste0('Not found any location match the sequence in ',chr.name[i]))
      }
    }
  }else{
    chr.seq <- getSeq(BS.hg19, chrMatch)
    seqr <- reverseComplement(seq)
    mtp <- matchPattern(pattern = seq, chr.seq)
    mtpr <- matchPattern(pattern = seqr, chr.seq)
    if(length(mtp)!=0){
      ranges <- mtp@ranges
      strand <- '+'
      location <- GRanges(seqnames = chrMatch,
                          ranges = ranges,
                          strand = strand)
      location.c <- changeIRanges(location, upstream, width)
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
                                 (location.c@ranges@width + location.c@ranges@start -1 ))
      print(paste0('The query sequence is sense strand according to the Genome ',chrMatch))
    }else if(length(mtpr)!=0){
      rangesr <- mtpr@ranges
      strand <- '-'
      location <- GRanges(seqnames = chrMatch,
                          ranges = rangesr,
                          strand = strand)
      
      location.c <- changeIRanges(location, upstream, width)
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
                                 (location.c@ranges@width + location.c@ranges@start ))
      print(paste0('The query sequence is anti-sense strand according to the BS.hg19 Genome ',chrMatch))
    }else{
      print(paste0('Not found any location match the sequence in ',chrMatch))
    }
  }
  

  if (!is.na(Flank.seq)) {
    if (!is.na(title)) {
      writeXStringSet(Flank.seq, filepath = paste0(title,'.fasta'))
    }
    return(Flank.seq)
  }else{
    return(NA)
  }

}



s.seq <- "GGAGGAAAAACTGTTTCATACAGAAGGCGT"

getSeqFgenome(s.seq)
