library(msa)
msaprintPDF <- function(align.obj = align.obj, title = title){
  cd <- getwd()
  msatemp <- file.path(substr(cd,1,2),'msatemp')
  
  dir.create(msatemp)
  outdir  <-  file.path(msatemp,"alignements")
  fasta_file  <-  paste(outdir,".fasta",sep="")
  
  output  <-  paste(outdir,"/",title,".tex", sep="")
  dir.create(outdir,showWarnings = FALSE)
  msaPrettyPrint(align.obj,
                 output="tex", 
                 showNames="left",
                 consensusColors = 'ColdHot',
                 askForOverwrite=FALSE, 
                 verbose=FALSE,
                 file = output, 
                 alFile = fasta_file)
  tools::texi2pdf(output, clean=TRUE)
  unlink(msatemp,recursive = T)
  file.copy(paste(title,'.pdf', sep = ""), paste('../',title,'.pdf', sep = ""))
  file.remove(paste(title,'.pdf', sep = ""))
}
