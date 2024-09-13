if (!require("BiocManager"))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("ggtree")

BiocManager::install("broadSeq")

abrowseVignettes(package = "broadSeq")

install.packages('devtools')

devtools::install_github("dasroy/broadSeq", build_vignettes = TRUE)
devtools::install_github("dasroy/broadSeq")

browseVignettes(package = "broadSeq")

BiocManager::install('BiocStyle')


install.packages('ggplot2')
install.packages('ggpubr')
