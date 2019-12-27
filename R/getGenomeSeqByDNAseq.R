library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
BS.hg19 <- BSgenome.Hsapiens.UCSC.hg19
changeIRanges <- function(granges.obj = granges.obj, 
                          upstream = integer, 
                          width = integer, ...)
{
  if (as.character(granges.obj@strand) == "-")
  {
    GRanges(seqnames = granges.obj@seqnames, 
            IRanges(end = (end(granges.obj@ranges) + upstream), 
                    start = (end(granges.obj@ranges) + upstream - width + 1)), 
            strand = granges.obj@strand, 
            seqinfo = granges.obj@seqinfo, 
            mcols(granges.obj))
  } else
  {
    GRanges(seqnames = granges.obj@seqnames, 
            IRanges(start = start(granges.obj@ranges) - upstream, 
                    width = width), strand = granges.obj@strand, 
            seqinfo = granges.obj@seqinfo, 
            mcols(granges.obj))
  }
}
getFlankGranges <- function(matchPattern.obj = matchPattern.obj, 
                            strand = strand, 
                            chromosome = chr, 
                            upstream = upstream, 
                            width = width, 
                            title = title, ...)
{
  ranges <- matchPattern.obj@ranges
  strand <- strand
  location <- GRanges(seqnames = chromosome, 
                      ranges = ranges, 
                      strand = strand)
  location.c <- changeIRanges(location, upstream, width)
  Flank_seq <- getSeq(BS.hg19, location.c)
  chrname <- as.character(location.c@seqnames)
  strandname <- as.character(location.c@strand)
  startname <- as.character(location.c@ranges@start)
  endname <- as.character((location.c@ranges@width + location.c@ranges@start - 1))
  names(Flank_seq) <- paste0(title, 
                             " ", 
                             chrname, 
                             " ", 
                             strandname, 
                             ":", 
                             " ", 
                             startname,  
                             "-", 
                             endname)
  if (strand == "+")
  {
    str <- "sense"
    print(paste0("The query sequence is ", 
                 str, 
                 " strand according to the BS.hg19 Genome ", 
                 chromosome))
  } else
  {
    str <- "anti-sense"
    print(paste0("The query sequence is ", 
                 str, 
                 " strand according to the BS.hg19 Genome ", 
                 chromosome))
  }
  return(Flank_seq)
}
getSeqFgenome <- function(seqstr = sequence, 
                          title = NA, 
                          chrMatch = NA, 
                          upstream = 0, 
                          width = NA, ...)
{
  seqstr <- seqstr
  if (is.na(width))
  {
    width <- nchar(seqstr)
  }
  seq <- DNAString(seqstr)
  Flank.seq <- NA
  if (is.na(chrMatch))
  {
    chr.name <- seqnames(BS.hg19)
    chr.name <- sample(chr.name, length(chr.name), replace = F)
    for (i in seq_len(length(chr.name)))
    {
      chr.seq <- getSeq(BS.hg19, chr.name[i])
      seqr <- reverseComplement(seq)
      mtp <- matchPattern(pattern = seq, chr.seq)
      mtpr <- matchPattern(pattern = seqr, chr.seq)
      if (length(mtp) != 0)
      {
        Flank.seq <- getFlankGranges(matchPattern.obj = mtp, 
                                     strand = "+", 
                                     chromosome = chr.name[i], 
                                     upstream = upstream, 
                                     width = width, 
                                     title = title)
        # break()
      } else if (length(mtpr) != 0)
      {
        Flank.seq <- getFlankGranges(matchPattern.obj = mtpr, 
                                     strand = "-", 
                                     chromosome = chr.name[i], 
                                     upstream = upstream, 
                                     width = width, 
                                     title = title)
        # break()
      } else
      {
        print(paste0("Not found any location match the sequence in ", chr.name[i]))
      }
    }
  } else
  {
    chr.seq <- getSeq(BS.hg19, chrMatch)
    seqr <- reverseComplement(seq)
    mtp <- matchPattern(pattern = seq, chr.seq)
    mtpr <- matchPattern(pattern = seqr, chr.seq)
    if (length(mtp) != 0)
    {
      Flank.seq <- getFlankGranges(matchPattern.obj = mtp, 
                                   strand = "+", 
                                   chromosome = chrMatch, 
                                   upstream = upstream, 
                                   width = width, 
                                   title = title)
      # break()
    } else if (length(mtpr) != 0)
    {
      Flank.seq <- getFlankGranges(matchPattern.obj = mtpr, 
                                   strand = "-", chromosome = chrMatch, 
                                   upstream = upstream, 
                                   width = width, 
                                   title = title)
      # break()
    } else
    {
      print(paste0("Not found any location match the sequence in ", chr.name[i]))
    }
  }
  if (!is.na(Flank.seq))
  {
    if (!is.na(title))
    {
      writeXStringSet(Flank.seq, filepath = paste0(title, ".fasta"))
    }
    return(Flank.seq)
  } else
  {
    return(NA)
  }
}





