Examples from the epilepsy dataset
==================================

Load packages, functions and data
---------------------------------

``` r
library(PopSV)
library(dplyr)
library(magrittr)
library(ggplot2)
library(GenomicRanges)
library(Gviz)
options(ucscChromosomeNames = FALSE)
source("EpiPopSV-scripts.R")

## Load genomic annotations
load("../data/gencodev19-gviz-simp.RData")
gtf.simp %<>% mutate(chromosome = gsub("chr", "", chromosome))
load("../data/SVdatabase.RData")
load("../data/genomicFeatures.RData")
## geneId -> geneName
gg = gen.feat.l$gene %>% select(geneId, geneName) %>% unique
geneidToGenename = gg$geneName
names(geneidToGenename) = gg$geneId
## Exons
exons.grl = list(exon = subset(gen.feat.l$exon, geneType == "protein_coding"))
gene.gr = makeGRangesFromDataFrame(gen.feat.l$gene, keep.extra.columns = TRUE)
## ExAC lof intolerant
exac = read.csv("../data/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", 
    as.is = TRUE, sep = "\t")
gene.intolerant = subset(exac, pLI > 0.9)$gene
exons.grl$exon.pli = subset(gen.feat.l$exon, geneName %in% gene.intolerant)
## Epilepsy genes
epilepsy.genes = scan("../data/EpilepsyGenes.txt", "a", quiet = TRUE)
exons.grl$exon.epi = subset(gen.feat.l$exon, geneName %in% epilepsy.genes)
## Convert to GRanges
exons.grl = lapply(exons.grl, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
## Enhancer map
load("../data/enhancerMap.RData")
eqtl.epi = subset(eqtl, gene_id %in% unique(exons.grl$exon.epi$geneId))
eqtl.epi$gene = geneidToGenename[eqtl.epi$gene_id]
seqlevels(dnaseMap) = gsub("chr", "", seqlevels(dnaseMap))
dnaseMap.epi = subset(dnaseMap, gene %in% epilepsy.genes)

## CNVs
cnv.all = read.table("../data/cnvs-PopSV-Epilepsy-198affected-301controls-5kb.tsv.gz", 
    header = TRUE, as.is = TRUE, sep = "\t")
cnv.all = cnv.all %>% mutate(project = ifelse(project == "affected", "epilepsy", 
    "parent"))
cnv.all.gr = cnv.all %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
cnv.all$prop.db = cnv.all.gr %>% dbProp(., svs.gr)
cnv.all$exon.epi = overlapsAny(cnv.all.gr, exons.grl$exon.epi)
cnv.all$enhancer.epi = overlapsAny(cnv.all.gr, dnaseMap.epi) | overlapsAny(cnv.all.gr, 
    eqtl.epi)
```

Example - PopSV calls and read coverage
---------------------------------------

``` r
plotCNVgene <- function(cnv.gr, bc.f, gene.gtf, range.gr, flanks = 0, flanks2 = 20000, 
    genes = NULL, ref.samples = NULL) {
    if (is.data.frame(cnv.gr)) 
        cnv.gr = makeGRangesFromDataFrame(cnv.gr, keep.extra.columns = TRUE)
    cnv.gr = subsetByOverlaps(cnv.gr, resize(range.gr, width(range.gr) + 2 * 
        flanks, fix = "center"))
    cnv.gr$type = ifelse(cnv.gr$z < 0, "deletion", "duplication")
    posl = min(start(cnv.gr), start(range.gr))
    posu = max(end(cnv.gr), end(range.gr))
    start(range.gr) = posl - flanks2
    end(range.gr) = posu + flanks2
    bc.df = read.bedix(bc.f, range.gr)
    bc.gr = makeGRangesFromDataFrame(bc.df, keep.extra.columns = TRUE)
    bc.cnv = bc.gr
    mcols(bc.cnv) = mcols(bc.cnv)[, unique(cnv.gr$sample)]
    if (!is.null(ref.samples)) 
        mcols(bc.gr) = mcols(bc.gr)[, ref.samples]
    ylims = c(min(as.matrix(mcols(bc.cnv))), max(as.matrix(mcols(bc.cnv))))
    ylims = c(min(ylims[1], apply(as.matrix(mcols(bc.gr)), 1, quantile, probs = 0.2)), 
        max(ylims[2], apply(as.matrix(mcols(bc.gr)), 1, quantile, probs = 0.8)))
    ylims = ylims + c(-1, 1) * diff(ylims) * 0.1
    bcref.t = DataTrack(bc.gr, type = "boxplot", name = "Coverage", do.out = FALSE, 
        ylim = ylims)
    bccnv.t = DataTrack(bc.cnv, groups = colnames(mcols(bc.cnv)), type = c("p", 
        "l"), legend = TRUE, lwd = 2, ylim = ylims)
    displayPars(bccnv.t) <- list(alpha.title = 1, alpha = 0.5)
    bcrefcnv.t = OverlayTrack(trackList = list(bcref.t, bccnv.t), ylim = ylims)
    cnv.t = AnnotationTrack(cnv.gr, group = cnv.gr$sample, feature = cnv.gr$type, 
        deletion = "steelblue", duplication = "indianred", name = "CNV", legend = TRUE)
    gene.gtf.gr = makeGRangesFromDataFrame(gene.gtf)
    gene.gtf = gene.gtf[overlapsAny(gene.gtf.gr, range.gr), ]
    if (!is.null(genes)) 
        gene.gtf = subset(gene.gtf, symbol %in% genes)
    axis.t <- GenomeAxisTrack()
    gene.t <- GeneRegionTrack(gene.gtf, genome = "hg19", chromosome = gene.gtf$chromosome[1], 
        name = "Gene", transcriptAnnotation = "symbol")
    plotTracks(list(axis.t, gene.t, cnv.t, bcrefcnv.t), groupAnnotation = "group", 
        sizes = c(1, 1, 2, 10))
    return(list(axis.t, gene.t, cnv.t, bcrefcnv.t))
}

bc.f = "../data/epilepsy-ref-bc-norm.tsv.bgz"

plotCNVgene(subset(cnv.all, project == "epilepsy"), bc.f, gtf.simp, subset(gene.gr, 
    gene == "CHD2"), genes = "CHD2", flanks = 1e+05)
plotCNVgene(subset(cnv.all, prop.db < 0.01 & project == "epilepsy" & (exon.epi | 
    enhancer.epi)), bc.f, gtf.simp, subset(gene.gr, gene == "BRAF"), genes = "BRAF", 
    flanks = 60000)
```
