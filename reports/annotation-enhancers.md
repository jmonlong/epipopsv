Enhancer annotations
====================

``` r
library(ggplot2)
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(data.table)
library(knitr)
library(tidyr)
```

We download eQTLs and DNase-to-promoter information that links non-coding regions to the gene they are most likely associated with.

GTEx eQTLs
----------

We downloaded GTEx eQTLs (v6p) from the [GTEx Portal](http://gtexportal.org/home/datasets). It requires registration (quick and free). We downloaded `GTEx_Analysis_v6p_eQTL.tar` tar file with *eGene and significant snp-gene associations based on permutations*.

``` r
gtex.f = untar("GTEx_Analysis_v6p_eQTL.tar", list = TRUE)
gtex.f = grep("snpgene_pairs", gtex.f, value = TRUE)
eqtl = lapply(gtex.f, function(fn) {
    tissue = gsub(".*/(.*)_Analysis.*", "\\1", fn)
    untar("GTEx_Analysis_v6p_eQTL.tar", fn)
    eqtl = fread(paste0("gunzip -c ", fn))
    file.remove(fn)
    eqtl = as.data.frame(eqtl)
    eqtl %>% mutate(chr = gsub("_.*", "", variant_id), start = gsub(".*_([0-9]+)_.*", 
        "\\1", variant_id), start = as.numeric(start), end = start + 1, score = -log10(pval_beta)) %>% 
        select(chr, start, end, gene_id, score) %>% mutate(tissue = tissue) %>% 
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
})
eqtl = do.call(c, eqtl)
eqtl %<>% as.data.frame %>% group_by(seqnames, start, end, gene_id) %>% summarize(score = mean(score), 
    nb.tis = n(), tissue = paste(tissue, collapse = ",")) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
```

DNase-to-promoter
-----------------

From Maurano et al ([Science 2012](http://www.sciencemag.org/content/337/6099/1190.full)), a distal hypersensitive site is linked to a promoter site if their presence are correlated across different cell types. All the correlated pairs of DNase sites are available [here](http://www.uwencode.org/proj/Science_Maurano_Humbert_et_al/data/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_35celltypeCategories.bed8.gz).

``` r
if (!file.exists("genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_35celltypeCategories.bed8.gz")) download.file("http://www.uwencode.org/proj/Science_Maurano_Humbert_et_al/data/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_35celltypeCategories.bed8.gz", 
    "genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_35celltypeCategories.bed8.gz")
dnaseMap.all = read.table("genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_35celltypeCategories.bed8.gz", 
    as.is = TRUE)
dnaseMap = dnaseMap.all[, 4:8]
colnames(dnaseMap) = c("gene", "chr", "start", "end", "score")
dnaseMap = makeGRangesFromDataFrame(dnaseMap, keep.extra.columns = TRUE)
```

Save the R objects
------------------

``` r
save(eqtl, dnaseMap, file = "../data/enhancerMap.RData")
```

We save these two types of information into a `.RData` file.
