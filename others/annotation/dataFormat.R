library(dplyr)
library(GenomicRanges)
library(data.table)

## Raw files downloaded following 'downloads.sh' script
ideo.f = "cytoBandIdeo.txt.gz"
gap.f = "gap.txt.gz"
gene.f = "gencode.v19.annotation.gene.gtf.gz"
exon.f = "gencode.v19.annotation.exon.gtf.gz"

## Centromere, telomere, gap annotation
chr.band = read.table(ideo.f, sep="\t", as.is=TRUE)
colnames(chr.band) = c("chr","start","end", "band", "type")
ct = read.table(gap.f,sep="\t", as.is=TRUE)
ct = ct[,c(2:4,8)]
colnames(ct) = c("chr","start","end","type")
## Adding telomeres as 10kbp from each chromosomal end
ct = rbind(ct, chr.band %>% group_by(chr) %>% summarize(start=min(start),end=1e4) %>% mutate(type="telomere"))
ct = rbind(ct, chr.band %>% group_by(chr) %>% summarize(start=max(end)-1e4,end=max(end)) %>% mutate(type="telomere"))
##
ct$chr = gsub("chr","",ct$chr)
centel.gr = with(ct,GRanges(chr,IRanges(start,end),type=type))
save(centel.gr, file="centelgap.RData")

## To parse GTF attributes
getAtt <- function(attributes, att.name="gene_id"){
  sub(paste0(".*",att.name," ([^;]+);.*"), "\\1", attributes)
}

## GENCODE Genes
genes = read.table(gene.f,as.is=TRUE, sep="\t")
genes = genes[,c(1,4,5,7,9)]
colnames(genes) = c("chr","start","end","strand","geneInfo")
genes$chr = gsub("chr","",genes$chr)
genes$geneId = getAtt(genes$geneInfo)
genes$geneType = getAtt(genes$geneInfo, "gene_type")
genes$geneName = getAtt(genes$geneInfo, "gene_name")
genes$geneInfo = NULL

## GENCODE Exons
exons = read.table(exon.f,as.is=TRUE, sep="\t")
exons = exons[,c(1,4,5,9)]
colnames(exons) = c("chr","start","end","geneInfo")
exons$chr = gsub("chr","",exons$chr)
exons$geneId = getAtt(exons$geneInfo)
exons$geneType = getAtt(exons$geneInfo, "gene_type")
exons$geneName = getAtt(exons$geneInfo, "gene_name")
exons$geneInfo = NULL

## Merge all annotations
gen.feat.l = list(gene=genes, exon=exons, centel=ct)
gen.feat.l = lapply(gen.feat.l, function(df) subset(df, chr %in% 1:22))
save(gen.feat.l, file="genomicFeatures.RData")
