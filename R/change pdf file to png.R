library(pdftools)
library(stringr)
filenames <- dir('../') %>% str_subset('[0-9].*\\.pdf')
for (i in seq_len(length(filenames))) {
  pdf <- paste0('../',filenames[i])
  png <- str_replace(paste0('../',filenames[i]), 'pdf','png')
  pdf_convert(pdf,
              dpi = 300,
              filenames = png)
}

