library(sangerseqR)

score.sanger <- function(sangerobj){
  traceMatrix <- tibble::as_tibble(traceMatrix(sangerobj))
  colnames(traceMatrix) <- c('A','C','G','T')
  traceMatrix
  data <- dplyr::mutate(traceMatrix, score = (traceMatrix$A < 50)+
                          (traceMatrix$C<50)+
                          (traceMatrix$G <50)+ 
                          (traceMatrix$`T`<50))
  data
}






