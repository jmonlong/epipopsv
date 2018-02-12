## Install QDNAseq
## source("http://bioconductor.org/biocLite.R")
## biocLite("QDNAseq")

library(QDNAseq)
library(BatchJobs) # to send jobs to HPC

bins <- getBinAnnotations(binSize=5)
save(bins, file='QDNAseq-bins.RData')

load('files.RData') # 'files.df' with a 'bam' column with the path to the bam file

## Read counts: one job per BAM file
reg <- makeRegistry(id="bcQDNAseq")
bcQDNAseq.f <- function(samp.ii, files.f, bins.f){
  library(QDNAseq)
  load(files.f)
  load(bins.f)
  binReadCounts(bins, bamfiles=files.df$bam[samp.ii], chunkSize=1e5, isProperPair=TRUE, minMapq=30)
}
batchMap(reg, bcQDNAseq.f, 1:nrow(files.df), more.args=list(files.f='files.RData', bins.f='QDNAseq-bins.RData'))
submitJobs(reg, findNotDone(reg), resources=list(walltime="30:00:00", nodes="1", cores="12", supervisor.group="bws-221-ae",queue='metaq'))
showStatus(reg)

## Correct read counts: one job per sample
reg2 <- makeRegistry(id="correctQDNAseq")
correctQDNAseq.f <- function(samp.ii){
  library(QDNAseq)
  library(BatchJobs)
  reg <- makeRegistry(id="bcQDNAseq")
  readCounts = loadResult(reg, samp.ii)
  readCountsFiltered <- applyFilters(readCounts)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  exportBins(copyNumbersSmooth, file=paste0('bc-corrected-',samp.ii,'.tsv'))
  return('Done')
}
batchMap(reg2, correctQDNAseq.f, findJobs(reg))
submitJobs(reg2, intersect(findDone(reg), findNotDone(reg2)), resources=list(walltime="1:00:00", nodes="1", cores="12", supervisor.group="bws-221-ae",queue='metaq'))
showStatus(reg2)

## Merge samples: one job
reg3 <- makeRegistry(id="merge")
merge.f <- function(outfile, samps){
  bc = lapply(samps, function(samp.ii){
    bc = read.table(paste0('bc-corrected-',samp.ii,'.tsv.gz'), as.is=TRUE, header=TRUE)
    bc = bc[,5, drop=FALSE]
    colnames(bc) = paste0('s',samp.ii)
    bc
  })
  bc.coords = read.table(paste0('bc-corrected-',samps[1],'.tsv'), as.is=TRUE, header=TRUE)
  bc.coords = bc.coords[,2:4]
  bc = cbind(bc.coords, do.call(cbind, bc))
  write.table(bc, file=outfile, row.names=FALSE, quote=FALSE, sep='\t')
  return('Done')
}
batchMap(reg3, merge.f, 'bc-corrected-QDNAseq-twins.tsv', more.args=list(samps=findDone(reg2)))
submitJobs(reg3, 1, resources=list(walltime="6:00:00", nodes="1", cores="12", supervisor.group="bws-221-ae",queue='metaq'))
showStatus(reg3)
