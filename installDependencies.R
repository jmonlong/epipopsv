pkgs = c("devtools","knitr","rmarkdown","magrittr","dplyr", "RColorBrewer","ggplot2","parallel", "tidyr", "ggdendro", "fpc", "formatR", "data.table")
pkgs.bioc = c("GenomicRanges","Gviz","BSgenome.Hsapiens.UCSC.hg19", "Rsamtools", "DNAcopy", "QDNAseq")
pkgs.github = c("jmonlong/PopSV")

install.packages(pkgs, repos="https://cloud.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=pkgs.bioc, suppressUpdates=TRUE)

sapply(pkgs.github, devtools::install_github)
