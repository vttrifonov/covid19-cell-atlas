#install.packages('renv')

#renv::init(bare=TRUE)
#renv::install('renv')

renv::install(c(
  'DT',
  'shiny',
  'data.table',
  'magrittr',
  'markdown',
  'rmarkdown',
  'Matrix',
  'bioc::limma',
  'reticulate',
  'languageserver',
  'httpgd',
  'ggplot2',
  'glue',
  'reshape2',
  'gridExtra',
  'bioc::fgsea',
  'kableExtra',
  'cowplot'
#  'readxl',
#  'mime'
))
