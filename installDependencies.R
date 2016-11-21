pkgs = c("devtools","knitr","magrittr","dplyr", "RColorBrewer","ggplot2","parallel", "tidyr", "ggdendro", "fpc")
pkgs.bioc = c("GenomicRanges","Gviz","BSgenome.Hsapiens.UCSC.hg19", "Rsamtools", "DNAcopy")
pgks.github = c("jmonlong/PopSV")

install.packages(pkgs, repos="https://cloud.r-project.org")

source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=pkgs.bioc, suppressUpdates=TRUE, suppressAutoUpdates=TRUE)

library(devtools)
sapply(pkgs.github, install_github(pkgs.github))
