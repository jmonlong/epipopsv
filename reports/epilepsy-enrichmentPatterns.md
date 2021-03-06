CNV enrichment patterns in epilepsy
===================================

Load packages, functions and data
---------------------------------

``` r
library(knitr)
library(RColorBrewer)
library(PopSV)
library(dplyr)
library(magrittr)
library(ggplot2)
library(GenomicRanges)
library(parallel)
library(tidyr)
source("EpiPopSV-scripts.R")

## Load genomic annotations
load("../data/SVdatabase.RData")
load("../data/genomicFeatures.RData")
## geneId -> geneName
gg = gen.feat.l$gene %>% select(geneId, geneName) %>% unique
geneidToGenename = gg$geneName
names(geneidToGenename) = gg$geneId
## Exons
exons.grl = list(exon = subset(gen.feat.l$exon, geneType == "protein_coding"))
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
## Load and annotate CNVs
cnv.all = read.table("../data/cnvs-PopSV-Epilepsy-198affected-301controls-5kb.tsv.gz", 
    header = TRUE, as.is = TRUE, sep = "\t")
cnv.all = cnv.all %>% mutate(project = ifelse(project == "affected", "patients", 
    "controls"))
cnv.all.gr = makeGRangesFromDataFrame(cnv.all, keep.extra.columns = TRUE)
cnv.all$prop.db = dbProp(cnv.all.gr, svs.gr)
cnv.all$prop.db.50 = cnv.all.gr %>% dbProp(., svs.gr, min.db.span = 0.5)  # 50% overlap
cnv.all$exon = overlapsAny(cnv.all.gr, exons.grl$exon)
epi.gr = reduce(subset(cnv.all.gr, project == "patients"))
cnv.all$exon.epi = overlapsAny(cnv.all.gr, exons.grl$exon.epi)
cnv.all$enhancer.epi = overlapsAny(cnv.all.gr, dnaseMap.epi) | overlapsAny(cnv.all.gr, 
    eqtl.epi)
epid = distanceToNearest(cnv.all.gr, exons.grl$exon.epi) %>% as.data.frame
cnv.all$exon.epi.d = Inf
cnv.all$exon.epi.d[epid$queryHits] = epid$distance
cnv.all$exon.epi.closest = NA
cnv.all$exon.epi.closest[epid$queryHits] = exons.grl$exon.epi$geneName[epid$subjectHits]
cnv.all %<>% group_by(project) %>% do(freq.range(., annotate.only = TRUE)) %>% 
    ungroup
info.df = cnv.all %>% select(sample, project) %>% unique
```

Also let's choose how many cores we want to use:

``` r
NB.CORES = 3
```

Exonic enrichment
-----------------

This part was run on a high-performance computing cluster. The code can be found in `epilepsy-enrichmentPatterns-HPC.R`.

Briefly, 150 samples were sub-sampled 100 times in the patient and in the control cohorts. For each sub-sampling, the CNVs are merged per cohort and the catalogs are overlapped with exonic sequences. All genes and genes with a predicted loss-of-function intolerance (pLi\>0.9) were used.

To control for the CNV size distribution, control regions were randomly selected in the reference genome controlling for the size and overlap with assembly gaps. To test for differences between patients and controls, CNVs from the two cohorts were permuted 10,000 times and we saved the difference between the medians of the two groups (see boxplots below).

We performed the analysis separately on *large* CNVs (\>=50 Kbp) and *small* CNVs (\<50 Kbp). We also run the analysis using rare CNVs (\<1% frequency in the external SV databases) only.

``` r
load("../data/permExCnv.RData")

enr.obs = lapply(permExCnv, function(l) l$enr.obs)
enr.obs = do.call(rbind, enr.obs)

enr.obs %<>% mutate(freq = factor(freq, levels = c("all", "rare"), labels = c("all CNVs", 
    "rare CNVs")), size = factor(size, levels = c("large", "small"), labels = c("> 50 Kbp", 
    "< 50 Kbp")))

ggplot(enr.obs, aes(x = feat, y = enr, fill = project)) + geom_boxplot(notch = TRUE) + 
    theme_bw() + facet_grid(size ~ freq, scales = "free") + xlab("") + ylab("fold-enrichment") + 
    scale_x_discrete(breaks = c("exon", "exon.pli"), labels = c("all genes", 
        "LoF\nintolerant\ngenes")) + geom_hline(yintercept = 1, linetype = 2) + 
    scale_fill_brewer(name = "", palette = "Set1") + theme(legend.position = c(0.01, 
    0.01), legend.justification = c(0, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
enr.exp.diff = lapply(permExCnv, function(l) l$enr.exp.diff)
enr.exp.diff = do.call(rbind, enr.exp.diff)
enr.exp.diff %<>% mutate(set = factor(paste(freq, size), levels = c("all large", 
    "rare large", "all small", "rare small"), labels = c("all CNVs > 50 Kbp", 
    "rare CNVs > 50 Kbp", "all CNVs < 50 Kbp", "rare CNVs < 50 Kbp")))

enr.obs.diff = enr.obs %>% mutate(set = factor(paste(freq, size), levels = c("all CNVs > 50 Kbp", 
    "rare CNVs > 50 Kbp", "all CNVs < 50 Kbp", "rare CNVs < 50 Kbp"))) %>% group_by(set, 
    feat, project) %>% summarize(enr = median(enr)) %>% group_by(set, feat) %>% 
    summarize(diff = enr[project == "patients"] - enr[project == "controls"])

ggplot(enr.exp.diff, aes(x = feat, y = diff)) + geom_hline(yintercept = 0, linetype = 2) + 
    geom_violin(fill = "grey90") + geom_point(data = enr.obs.diff, size = 5, 
    color = "red") + theme_bw() + facet_wrap(~set, scales = "free", ncol = 2) + 
    xlab("") + ylab("difference between patients and controls") + scale_x_discrete(breaks = c("exon", 
    "exon.pli"), labels = c("all genes", "LoF\nintolerant\ngenes")) + scale_fill_brewer(name = "", 
    palette = "Set1")
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-3-2.png)

*The grey violin plot shows the distribution of the differences across the 10,000 permutations. The red dot represents the difference observed between the patients and controls.*

The empirical P-value from the 10,000 permutations:

``` r
enr.obs.diff %>% dplyr::rename(diff.obs = diff) %>% merge(enr.exp.diff) %>% 
    group_by(set, feat) %>% summarize(pv.enr = (sum(diff >= diff.obs) + 1)/(n() + 
    1), pv.dep = (sum(diff <= diff.obs) + 1)/(n() + 1)) %>% kable
```

| set                 | feat     |     pv.enr|     pv.dep|
|:--------------------|:---------|----------:|----------:|
| all CNVs \> 50 Kbp  | exon     |  0.2154785|  0.7846215|
| all CNVs \> 50 Kbp  | exon.pli |  0.0001000|  1.0000000|
| rare CNVs \> 50 Kbp | exon     |  0.0001000|  1.0000000|
| rare CNVs \> 50 Kbp | exon.pli |  0.0001000|  1.0000000|
| all CNVs \< 50 Kbp  | exon     |  0.6210379|  0.3790621|
| all CNVs \< 50 Kbp  | exon.pli |  0.0001000|  1.0000000|
| rare CNVs \< 50 Kbp | exon     |  0.0001000|  1.0000000|
| rare CNVs \< 50 Kbp | exon.pli |  0.0001000|  1.0000000|

Rare exonic CNVs are less private in the epilepsy cohort
--------------------------------------------------------

### All epilepsy patients, controls down-sampled

``` r
rare.ex = cnv.all %>% filter(prop.db < 0.01, exon)
rare.ex.f = rare.ex %>% group_by(project) %>% do(freqSS(.))
rare.ex.f %>% group_by(project, nb, rep) %>% summarize(n = n()/nb[1]) %>% group_by(project, 
    rep) %>% arrange(desc(nb)) %>% mutate(cn = cumsum(n), cprop = cn/sum(n)) %>% 
    filter(nb > 1, nb < 11) %>% group_by(project, nb) %>% summarize(cprop.5 = quantile(cprop, 
    probs = 0.05), cprop.95 = quantile(cprop, probs = 0.95), cprop = median(cprop)) %>% 
    ggplot(aes(x = nb, y = cprop, colour = project, fill = project)) + geom_ribbon(aes(ymin = cprop.5, 
    ymax = cprop.95), colour = FALSE, alpha = 0.2) + geom_line(alpha = 0.7) + 
    geom_point(size = 2) + scale_colour_brewer(name = "", palette = "Set1") + 
    scale_fill_brewer(name = "", palette = "Set1") + theme_bw() + scale_x_continuous(breaks = 2:10, 
    labels = paste0(c(2:10), "+")) + xlab("CNV recurrence") + ylab("proportion of rare exonic CNVs") + 
    theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-5-1.png)

### Same after exteme samples removal

In each cohort, let's remove the top 20 samples with the most rare exonic non-private CNVs.

``` r
rare.epi.samp.o = rare.ex.f %>% filter(nb > 1) %>% group_by(project, sample) %>% 
    summarize(cnv = n()) %>% arrange(desc(cnv)) %>% group_by(project) %>% do(head(., 
    20))
rare.ex.f.noext = rare.ex %>% filter(!(sample %in% rare.epi.samp.o$sample)) %>% 
    group_by(project) %>% do(freqSS(., sample.size = 178))
rare.ex.f.noext %>% group_by(project, nb, rep) %>% summarize(n = n()/nb[1]) %>% 
    group_by(project, rep) %>% arrange(desc(nb)) %>% mutate(cn = cumsum(n), 
    cprop = cn/sum(n)) %>% filter(nb > 1, nb < 11) %>% group_by(project, nb) %>% 
    summarize(cprop.5 = quantile(cprop, probs = 0.05), cprop.95 = quantile(cprop, 
        probs = 0.95), cprop = median(cprop)) %>% ggplot(aes(x = nb, y = cprop, 
    colour = project, fill = project)) + geom_ribbon(aes(ymin = cprop.5, ymax = cprop.95), 
    colour = FALSE, alpha = 0.2) + geom_line(alpha = 0.7) + geom_point(size = 2) + 
    scale_colour_brewer(name = "", palette = "Set1") + scale_fill_brewer(name = "", 
    palette = "Set1") + theme_bw() + scale_x_continuous(breaks = 2:10, labels = paste0(c(2:10), 
    "+")) + xlab("CNV recurrence") + ylab("proportion of rare exonic CNVs") + 
    ggtitle("Top 20 most extreme samples removed") + theme(legend.position = c(0.99, 
    0.99), legend.justification = c(1, 1))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-6-1.png)

### French-Canadians only

``` r
french.canadians = scan("../data/frenchCanadiansControls.txt", "a", quiet = TRUE)
cnv.fc = cnv.all %>% filter(project == "patients" | sample %in% french.canadians)
fc.ss = cnv.fc %>% select(sample, project) %>% unique %>% group_by(project) %>% 
    summarize(n = n()) %>% .$n %>% unlist %>% min
rare.ex.fc = cnv.fc %>% filter(exon, prop.db < 0.01) %>% group_by(project) %>% 
    do(freqSS(., sample.size = fc.ss))
rare.ex.fc %>% group_by(project, nb, rep) %>% summarize(n = n()/nb[1]) %>% group_by(project, 
    rep) %>% arrange(desc(nb)) %>% mutate(cn = cumsum(n), cprop = cn/sum(n)) %>% 
    filter(nb > 1, nb < 11) %>% group_by(project, nb) %>% summarize(cprop.5 = quantile(cprop, 
    probs = 0.05), cprop.95 = quantile(cprop, probs = 0.95), cprop = median(cprop)) %>% 
    ggplot(aes(x = nb, y = cprop, colour = project, fill = project)) + geom_ribbon(aes(ymin = cprop.5, 
    ymax = cprop.95), colour = FALSE, alpha = 0.2) + geom_line(alpha = 0.7) + 
    geom_point(size = 2) + scale_colour_brewer(name = "", palette = "Set1") + 
    scale_fill_brewer(name = "", palette = "Set1") + theme_bw() + scale_x_continuous(breaks = 2:10, 
    labels = paste0(c(2:10), "+")) + xlab("CNV recurrence") + ylab("proportion of rare exonic CNVs") + 
    ggtitle("French-Canadians only") + theme(legend.position = c(0.99, 0.99), 
    legend.justification = c(1, 1))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-7-1.png)

### Permutation test

``` r
diff.obs = rare.ex.f %>% group_by(project, nb) %>% summarize(n = n()/nb[1]) %>% 
    group_by(project) %>% summarize(non.priv = sum(n[nb > 1])/sum(n)) %>% arrange(project) %>% 
    .$non.priv %>% diff
diff.exp = replicate(1000, {
    rare.ex.f %>% ungroup %>% mutate(project = sample(project)) %>% group_by(project, 
        nb) %>% summarize(n = n()/nb[1]) %>% group_by(project) %>% summarize(non.priv = sum(n[nb > 
        1])/sum(n)) %>% arrange(project) %>% .$non.priv %>% diff
})

diff.obs.noext = rare.ex.f.noext %>% group_by(project, nb) %>% summarize(n = n()/nb[1]) %>% 
    group_by(project) %>% summarize(non.priv = sum(n[nb > 1])/sum(n)) %>% arrange(project) %>% 
    .$non.priv %>% diff
diff.exp.noext = replicate(1000, {
    rare.ex.f.noext %>% ungroup %>% mutate(project = sample(project)) %>% group_by(project, 
        nb) %>% summarize(n = n()/nb[1]) %>% group_by(project) %>% summarize(non.priv = sum(n[nb > 
        1])/sum(n)) %>% arrange(project) %>% .$non.priv %>% diff
})

diff.obs.fc = rare.ex.fc %>% group_by(project, nb) %>% summarize(n = n()/nb[1]) %>% 
    group_by(project) %>% summarize(non.priv = sum(n[nb > 1])/sum(n)) %>% arrange(project) %>% 
    .$non.priv %>% diff
diff.exp.fc = replicate(1000, {
    rare.ex.fc %>% ungroup %>% mutate(project = sample(project)) %>% group_by(project, 
        nb) %>% summarize(n = n()/nb[1]) %>% group_by(project) %>% summarize(non.priv = sum(n[nb > 
        1])/sum(n)) %>% arrange(project) %>% .$non.priv %>% diff
})
```

``` r
priv.test = data.frame(test = c("all samples", "top 20 extreme samples removed", 
    "french-canadian"), pv = c((1 + sum(diff.obs <= diff.exp))/(1 + length(diff.exp)), 
    (1 + sum(diff.obs.noext <= diff.exp.noext))/(1 + length(diff.exp.noext)), 
    (1 + sum(diff.obs.fc <= diff.exp.fc))/(1 + length(diff.exp.fc))))
kable(priv.test)
```

| test                           |        pv|
|:-------------------------------|---------:|
| all samples                    |  0.000999|
| top 20 extreme samples removed |  0.000999|
| french-canadian                |  0.001998|

Non-coding rare CNVs
--------------------

``` r
freqEpiRareSS <- function(df, sample.size = 198, nb.rep = 100) {
    if (length(unique(df$sample)) > sample.size) {
        res = mclapply(1:nb.rep, function(ii) {
            df %>% filter(sample %in% sample(unique(sample), sample.size)) %>% 
                filter(exon.epi.d < 3e+06) %>% do(freq.range(., annotate.only = TRUE)) %>% 
                filter(exon.epi.d < 4e+05, prop.db < 0.01) %>% mutate(rep = ii) %>% 
                as.data.frame
        }, mc.cores = NB.CORES)
        return(do.call(rbind, res))
    } else {
        return(df %>% freq.range(annotate.only = TRUE) %>% filter(exon.epi.d < 
            4e+05, prop.db < 0.01) %>% mutate(rep = 1))
    }
}
cnv.ss = cnv.all %>% group_by(project) %>% do(freqEpiRareSS(., nb.rep = 100))
```

### All rare non-coding CNVs

``` r
cnb.nc = cnv.ss %>% filter(prop.db < 0.01, nb < 10, !exon) %>% group_by(project, 
    rep) %>% arrange(exon.epi.d) %>% filter(!duplicated(sample)) %>% filter(exon.epi.d < 
    3e+05)
cnb.nc %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% group_by(project, 
    d) %>% summarize(cnb.5 = quantile(cnb, probs = 0.05), cnb.95 = quantile(cnb, 
    probs = 0.95), cnb = median(cnb)) %>% filter(d < 3e+05) %>% ggplot(aes(x = d/1000, 
    y = cnb, colour = project, fill = project)) + geom_line(size = 2) + geom_ribbon(aes(ymax = cnb.95, 
    ymin = cnb.5), linetype = 2, alpha = 0.2) + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("cumulative affected samples") + scale_fill_brewer(name = "", palette = "Set1") + 
    scale_colour_brewer(name = "", palette = "Set1") + guides(alpha = FALSE) + 
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
cnb.nc.or = cnb.nc %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% 
    group_by(project, d) %>% summarize(cnb = median(cnb)) %>% filter(d < 3e+05) %>% 
    spread(project, cnb) %>% mutate(odds.ratio = (patients/(198 - patients))/(controls/(198 - 
    controls)))
ggplot(cnb.nc.or, aes(x = d/1000, y = log(odds.ratio))) + geom_area(alpha = 0.1) + 
    geom_line() + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("log odds ratio") + geom_hline(yintercept = 0, linetype = 2)
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-11-2.png)

``` r
cnb.nc.or %>% filter(d %in% c(5000, 50000, 1e+05)) %>% kable
```

|      d|  controls|  patients|  odds.ratio|
|------:|---------:|---------:|-----------:|
|  5e+03|        17|        30|    1.901261|
|  5e+04|        62|        76|    1.366473|
|  1e+05|        88|       101|    1.301546|

``` r
ks.test(cnb.nc$exon.epi.d[which(cnb.nc$project == "controls")], cnb.nc$exon.epi.d[which(cnb.nc$project != 
    "controls")])
```

    ## 
    ##  Two-sample Kolmogorov-Smirnov test
    ## 
    ## data:  cnb.nc$exon.epi.d[which(cnb.nc$project == "controls")] and cnb.nc$exon.epi.d[which(cnb.nc$project != "controls")]
    ## D = 0.13944, p-value = 0.00533
    ## alternative hypothesis: two-sided

``` r
wilcox.test(cnb.nc$exon.epi.d[which(cnb.nc$project == "controls")], cnb.nc$exon.epi.d[which(cnb.nc$project != 
    "controls")])
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  cnb.nc$exon.epi.d[which(cnb.nc$project == "controls")] and cnb.nc$exon.epi.d[which(cnb.nc$project != "controls")]
    ## W = 1247900, p-value = 0.06973
    ## alternative hypothesis: true location shift is not equal to 0

#### Deletions

``` r
cnb.nc.del = cnv.ss %>% filter(z < 0, prop.db < 0.01, nb < 10, !exon) %>% group_by(project, 
    rep) %>% arrange(exon.epi.d) %>% filter(!duplicated(sample)) %>% filter(exon.epi.d < 
    3e+05)
cnb.nc.del %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% group_by(project, 
    d) %>% summarize(cnb.5 = quantile(cnb, probs = 0.05), cnb.95 = quantile(cnb, 
    probs = 0.95), cnb = median(cnb)) %>% filter(d < 3e+05) %>% ggplot(aes(x = d/1000, 
    y = cnb, colour = project, fill = project)) + geom_line(size = 2) + geom_ribbon(aes(ymax = cnb.95, 
    ymin = cnb.5), linetype = 2, alpha = 0.2) + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("cumulative affected samples") + scale_fill_brewer(name = "", palette = "Set1") + 
    scale_colour_brewer(name = "", palette = "Set1") + guides(alpha = FALSE) + 
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/ncdel-1.png)

``` r
cnb.nc.del.or = cnb.nc.del %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% 
    group_by(project, d) %>% summarize(cnb = median(cnb)) %>% filter(d < 3e+05) %>% 
    spread(project, cnb) %>% mutate(odds.ratio = (patients/(198 - patients))/(controls/(198 - 
    controls)))
ggplot(cnb.nc.del.or, aes(x = d/1000, y = log(odds.ratio))) + geom_area(alpha = 0.1) + 
    geom_line() + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("log odds ratio") + geom_hline(yintercept = 0, linetype = 2)
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/ncdel-2.png)

``` r
cnb.nc.del.or %>% filter(d %in% c(5000, 50000, 1e+05)) %>% kable
```

|      d|  controls|  patients|  odds.ratio|
|------:|---------:|---------:|-----------:|
|  5e+03|        10|        21|    2.230509|
|  5e+04|        38|        51|    1.460795|
|  1e+05|        52|        69|    1.501789|

``` r
ks.test(cnb.nc.del$exon.epi.d[which(cnb.nc.del$project == "controls")], cnb.nc.del$exon.epi.d[which(cnb.nc.del$project != 
    "controls")])
```

    ## 
    ##  Two-sample Kolmogorov-Smirnov test
    ## 
    ## data:  cnb.nc.del$exon.epi.d[which(cnb.nc.del$project == "controls")] and cnb.nc.del$exon.epi.d[which(cnb.nc.del$project != "controls")]
    ## D = 0.14261, p-value = 0.02305
    ## alternative hypothesis: two-sided

``` r
wilcox.test(cnb.nc.del$exon.epi.d[which(cnb.nc.del$project == "controls")], 
    cnb.nc.del$exon.epi.d[which(cnb.nc.del$project != "controls")])
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  cnb.nc.del$exon.epi.d[which(cnb.nc.del$project == "controls")] and cnb.nc.del$exon.epi.d[which(cnb.nc.del$project != "controls")]
    ## W = 573640, p-value = 0.1291
    ## alternative hypothesis: true location shift is not equal to 0

To go further we test each epilepsy gene to see if it has more rare non-coding deletions in the epilepsy cohort versus the controls.

``` r
cnv.all %>% filter(!exon, prop.db < 0.01, exon.epi.d < 2e+05, z < 0) %>% group_by(exon.epi.closest, 
    project) %>% summarize(sample = length(unique(sample))) %>% spread(project, 
    sample, fill = 0) %>% mutate(pv = fisher.test(cbind(c(301 - controls, controls), 
    c(198 - patients, patients)))$p.value) %>% arrange(pv) %>% ungroup %>% head(10) %>% 
    kable
```

| exon.epi.closest |  controls|  patients|         pv|
|:-----------------|---------:|---------:|----------:|
| GABRD            |         0|         4|  0.0243363|
| DNM1             |         0|         3|  0.0619015|
| NGLY1            |         0|         3|  0.0619015|
| SERPINI1         |         0|         3|  0.0619015|
| RYR2             |         1|         4|  0.0835300|
| GRIN2B           |         2|         5|  0.1198971|
| B3GNT4           |         0|         2|  0.1569645|
| CACNB4           |         0|         2|  0.1569645|
| SLC2A1           |         0|         2|  0.1569645|
| RBFOX1           |         9|        11|  0.1671548|

#### Duplications

``` r
cnb.nc.dup = cnv.ss %>% filter(z > 0, prop.db < 0.01, nb < 10, !exon) %>% group_by(project, 
    rep) %>% arrange(exon.epi.d) %>% filter(!duplicated(sample)) %>% filter(exon.epi.d < 
    3e+05)
cnb.nc.dup %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% group_by(project, 
    d) %>% summarize(cnb.5 = quantile(cnb, probs = 0.05), cnb.95 = quantile(cnb, 
    probs = 0.95), cnb = median(cnb)) %>% filter(d < 3e+05) %>% ggplot(aes(x = d/1000, 
    y = cnb, colour = project, fill = project)) + geom_line(size = 2) + geom_ribbon(aes(ymax = cnb.95, 
    ymin = cnb.5), linetype = 2, alpha = 0.2) + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("cumulative affected samples") + scale_fill_brewer(name = "", palette = "Set1") + 
    scale_colour_brewer(name = "", palette = "Set1") + guides(alpha = FALSE) + 
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/ncdup-1.png)

``` r
cnb.nc.dup.or = cnb.nc.dup %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% 
    group_by(project, d) %>% summarize(cnb = median(cnb)) %>% filter(d < 3e+05) %>% 
    spread(project, cnb) %>% mutate(odds.ratio = (patients/(198 - patients))/(controls/(198 - 
    controls)))
ggplot(cnb.nc.dup.or, aes(x = d/1000, y = log(odds.ratio))) + geom_area(alpha = 0.1) + 
    geom_line() + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("log odds ratio") + geom_hline(yintercept = 0, linetype = 2)
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/ncdup-2.png)

``` r
cnb.nc.dup.or %>% filter(d %in% c(5000, 50000, 1e+05)) %>% kable
```

|      d|  controls|  patients|  odds.ratio|
|------:|---------:|---------:|-----------:|
|  5e+03|         7|        10|    1.451368|
|  5e+04|        33|        35|    1.073620|
|  1e+05|        52|        54|    1.052885|

``` r
ks.test(cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project == "controls")], cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project != 
    "controls")])
```

    ## 
    ##  Two-sample Kolmogorov-Smirnov test
    ## 
    ## data:  cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project == "controls")] and cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project != "controls")]
    ## D = 0.11698, p-value = 0.1604
    ## alternative hypothesis: two-sided

``` r
wilcox.test(cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project == "controls")], 
    cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project != "controls")])
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project == "controls")] and cnb.nc.dup$exon.epi.d[which(cnb.nc.dup$project != "controls")]
    ## W = 506430, p-value = 0.2476
    ## alternative hypothesis: true location shift is not equal to 0

### Rare non-coding CNVs that overlaps functional annotations

Focusing on CNVs that overlaps an eQTL for the epilepsy gene, or a DNase I hypersensitive site associated to the promoter of the epilepsy gene.

``` r
cnb.df = cnv.ss %>% filter(prop.db < 0.01, nb < 10, !exon, enhancer.epi) %>% 
    group_by(project, rep) %>% arrange(exon.epi.d) %>% filter(!duplicated(sample)) %>% 
    filter(exon.epi.d < 3e+05)
cnb.df %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% group_by(project, 
    d) %>% summarize(cnb.5 = quantile(cnb, probs = 0.05), cnb.95 = quantile(cnb, 
    probs = 0.95), cnb = median(cnb)) %>% filter(d < 3e+05) %>% ggplot(aes(x = d/1000, 
    y = cnb, colour = project, fill = project)) + geom_line(size = 2) + geom_ribbon(aes(ymax = cnb.95, 
    ymin = cnb.5), alpha = 0.2, linetype = 2) + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("cumulative affected samples") + scale_fill_brewer(name = "", palette = "Set1") + 
    scale_colour_brewer(name = "", palette = "Set1") + theme(legend.position = c(0.99, 
    0.01), legend.justification = c(1, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
cnb.or = cnb.df %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% 
    group_by(project, d) %>% summarize(cnb = median(cnb)) %>% filter(d < 3e+05) %>% 
    spread(project, cnb) %>% mutate(odds.ratio = (patients/(198 - patients))/(controls/(198 - 
    controls)))
ggplot(cnb.or, aes(x = d/1000, y = log(odds.ratio))) + geom_area(alpha = 0.1) + 
    geom_line() + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("log odds ratio") + geom_hline(yintercept = 0, linetype = 2)
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-12-2.png)

``` r
cnb.or %>% filter(d %in% c(5000, 50000, 1e+05)) %>% kable
```

|      d|  controls|  patients|  odds.ratio|
|------:|---------:|---------:|-----------:|
|  5e+03|         8|        19|    2.520950|
|  5e+04|        31|        45|    1.584440|
|  1e+05|        49|        59|    1.290706|

``` r
ks.test(cnb.df$exon.epi.d[which(cnb.df$project == "controls")], cnb.df$exon.epi.d[which(cnb.df$project != 
    "controls")])
```

    ## 
    ##  Two-sample Kolmogorov-Smirnov test
    ## 
    ## data:  cnb.df$exon.epi.d[which(cnb.df$project == "controls")] and cnb.df$exon.epi.d[which(cnb.df$project != "controls")]
    ## D = 0.25554, p-value = 0.0001227
    ## alternative hypothesis: two-sided

``` r
wilcox.test(cnb.df$exon.epi.d[which(cnb.df$project == "controls")], cnb.df$exon.epi.d[which(cnb.df$project != 
    "controls")])
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  cnb.df$exon.epi.d[which(cnb.df$project == "controls")] and cnb.df$exon.epi.d[which(cnb.df$project != "controls")]
    ## W = 362060, p-value = 3.692e-05
    ## alternative hypothesis: true location shift is not equal to 0

``` r
cnb.rec = cnv.ss %>% filter(prop.db < 0.01, nb > 1, nb < 10, !exon, enhancer.epi) %>% 
    group_by(project, rep) %>% arrange(exon.epi.d) %>% filter(!duplicated(sample)) %>% 
    filter(exon.epi.d < 3e+05)
cnb.rec %>% do(decdf(., "exon.epi.d")) %>% group_by(project, d) %>% summarize(cnb.5 = quantile(cnb, 
    probs = 0.05), cnb.95 = quantile(cnb, probs = 0.95), cnb = median(cnb)) %>% 
    filter(d < 3e+05) %>% ggplot(aes(x = d/1000, y = cnb, colour = project, 
    fill = project)) + geom_line(size = 2) + geom_ribbon(aes(ymax = cnb.95, 
    ymin = cnb.5), alpha = 0.2, linetype = 2) + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("cumulative affected samples") + ggtitle("Non-private") + scale_fill_brewer(name = "", 
    palette = "Set1") + scale_colour_brewer(name = "", palette = "Set1") + theme(legend.position = c(0.99, 
    0.01), legend.justification = c(1, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-12-3.png)

``` r
cnb.rec.or = cnb.rec %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% 
    group_by(project, d) %>% summarize(cnb = median(cnb)) %>% filter(d < 3e+05) %>% 
    spread(project, cnb) %>% mutate(odds.ratio = (patients/(198 - patients))/(controls/(198 - 
    controls)))
ggplot(cnb.rec.or, aes(x = d/1000, y = log(odds.ratio))) + geom_area(alpha = 0.1) + 
    geom_line() + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("log odds ratio") + geom_hline(yintercept = 0, linetype = 2)
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-12-4.png)

``` r
cnb.rec.or %>% filter(d %in% c(5000, 50000, 1e+05)) %>% kable
```

|      d|  controls|  patients|  odds.ratio|
|------:|---------:|---------:|-----------:|
|  5e+03|         2|        10|    5.212766|
|  5e+04|         4|        17|    4.555249|
|  1e+05|         9|        19|    2.229050|

``` r
ks.test(cnb.rec$exon.epi.d[which(cnb.rec$project == "controls")], cnb.rec$exon.epi.d[which(cnb.rec$project != 
    "controls")])
```

    ## 
    ##  Two-sample Kolmogorov-Smirnov test
    ## 
    ## data:  cnb.rec$exon.epi.d[which(cnb.rec$project == "controls")] and cnb.rec$exon.epi.d[which(cnb.rec$project != "controls")]
    ## D = 0.60355, p-value = 5.373e-07
    ## alternative hypothesis: two-sided

``` r
wilcox.test(cnb.rec$exon.epi.d[which(cnb.rec$project == "controls")], cnb.rec$exon.epi.d[which(cnb.rec$project != 
    "controls")])
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  cnb.rec$exon.epi.d[which(cnb.rec$project == "controls")] and cnb.rec$exon.epi.d[which(cnb.rec$project != "controls")]
    ## W = 31706, p-value = 1.255e-06
    ## alternative hypothesis: true location shift is not equal to 0

#### Deletions

``` r
cnb.enh.del = cnv.ss %>% filter(z < 0, prop.db < 0.01, nb < 10, !exon, enhancer.epi) %>% 
    group_by(project, rep) %>% arrange(exon.epi.d) %>% filter(!duplicated(sample)) %>% 
    filter(exon.epi.d < 3e+05)
cnb.enh.del %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% group_by(project, 
    d) %>% summarize(cnb.5 = quantile(cnb, probs = 0.05), cnb.95 = quantile(cnb, 
    probs = 0.95), cnb = median(cnb)) %>% filter(d < 3e+05) %>% ggplot(aes(x = d/1000, 
    y = cnb, colour = project, fill = project)) + geom_line(size = 2) + geom_ribbon(aes(ymax = cnb.95, 
    ymin = cnb.5), linetype = 2, alpha = 0.2) + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("cumulative affected samples") + scale_fill_brewer(name = "", palette = "Set1") + 
    scale_colour_brewer(name = "", palette = "Set1") + guides(alpha = FALSE) + 
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/ncdelenh-1.png)

#### Duplications

``` r
cnb.enh.dup = cnv.ss %>% filter(z > 0, prop.db < 0.01, nb < 10, !exon, enhancer.epi) %>% 
    group_by(project, rep) %>% arrange(exon.epi.d) %>% filter(!duplicated(sample)) %>% 
    filter(exon.epi.d < 3e+05)
cnb.enh.dup %>% group_by(project, rep) %>% do(decdf(., "exon.epi.d")) %>% group_by(project, 
    d) %>% summarize(cnb.5 = quantile(cnb, probs = 0.05), cnb.95 = quantile(cnb, 
    probs = 0.95), cnb = median(cnb)) %>% filter(d < 3e+05) %>% ggplot(aes(x = d/1000, 
    y = cnb, colour = project, fill = project)) + geom_line(size = 2) + geom_ribbon(aes(ymax = cnb.95, 
    ymin = cnb.5), linetype = 2, alpha = 0.2) + theme_bw() + xlab("distance to nearest epilepsy exon (kb)") + 
    ylab("cumulative affected samples") + scale_fill_brewer(name = "", palette = "Set1") + 
    scale_colour_brewer(name = "", palette = "Set1") + guides(alpha = FALSE) + 
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/ncdupenh-1.png)

Rare CNVs hitting epilepsy genes
--------------------------------

### Controlling for the gene size

Epilepsy genes tend to be large, and larger genes are more likely hit by variants (including CNVs). We can select genes with similar gene size in order to control for the size effect. We check that in the distribution of gene size are satisfactory.

``` r
exon.cnv.genes = cnv.all %>% filter(project == "patients") %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    subsetByOverlaps(exons.grl$exon, .) %>% mcols %>% .$geneName %>% unique
gene.sum = exons.grl$exon %>% as.data.frame %>% group_by(geneName) %>% summarize(exon.size = sum(end - 
    start), start = min(start), end = max(end), size = end - start, nb.exon = n()) %>% 
    mutate(epilepsy = geneName %in% epilepsy.genes, cnv = geneName %in% exon.cnv.genes)
size.bk = c(quantile(subset(gene.sum, cnv)$size, probs = seq(0, 1, 0.2)), Inf)
gene.sum$size.class = cut(gene.sum$size, size.bk, include.lowest = TRUE)
gene.sum.cnv = subset(gene.sum, cnv)
gene.cont.ii = lapply(levels(gene.sum.cnv$size.class), function(sc) {
    sample(which(gene.sum$size.class == sc), sum(gene.sum.cnv$size.class == 
        sc))
})
gene.cont = gene.sum$geneName[unlist(gene.cont.ii)]
gene.sum.n = rbind(gene.sum %>% mutate(set = "all genes"), gene.sum %>% filter(geneName %in% 
    gene.cont) %>% mutate(set = "CNV-hit genes control"), gene.sum %>% filter(epilepsy) %>% 
    mutate(set = "epilepsy genes"), gene.sum %>% filter(cnv) %>% mutate(set = "CNV-hit genes"))

colpal = brewer.pal(3, "Set1")[c(1, 2, 2, 3)]
ggplot(gene.sum.n, aes(x = size, colour = set, linetype = set)) + stat_density(size = 2, 
    geom = "line", position = "dodge") + theme_bw() + scale_x_log10() + scale_linetype_manual(values = c(1, 
    1, 2, 1), name = "") + scale_colour_manual(values = colpal, name = "") + 
    ylab("density") + xlab("gene size (bp)") + theme(legend.position = c(0, 
    1), legend.justification = c(0, 1))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-13-1.png)

Eventially we can check that the exonic sequence size and number of exons are also controlled by this approach.

``` r
ggplot(gene.sum.n, aes(x = exon.size, colour = set, linetype = set)) + stat_density(size = 2, 
    geom = "line", position = "dodge") + theme_bw() + scale_x_log10() + scale_linetype_manual(values = c(1, 
    1, 2, 1), name = "") + scale_colour_manual(values = colpal, name = "") + 
    ylab("density") + xlab("exonic size (bp)") + theme(legend.position = c(0, 
    1), legend.justification = c(0, 1))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
ggplot(gene.sum.n, aes(x = nb.exon, colour = set, linetype = set)) + stat_density(size = 2, 
    geom = "line", position = "dodge") + theme_bw() + scale_x_log10() + scale_linetype_manual(values = c(1, 
    1, 2, 1), name = "") + scale_colour_manual(values = colpal, name = "") + 
    ylab("density") + xlab("number of exons") + theme(legend.position = c(0, 
    1), legend.justification = c(0, 1))
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-14-2.png)

Then the new sampling-based test:

``` r
testGenesCnvs <- function(cnv.df, nb.perm = 1000, nb.class = 5) {
    exon.cnv.genes = cnv.df %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
        subsetByOverlaps(exons.grl$exon, .) %>% mcols %>% .$geneName %>% unique
    gene.sum %<>% mutate(cnv = geneName %in% exon.cnv.genes)
    size.bk = c(quantile(subset(gene.sum, cnv)$size, probs = seq(0, 1, 1/nb.class)), 
        Inf)
    gene.sum$size.class = cut(gene.sum$size, size.bk, include.lowest = TRUE)
    gene.sum.cnv = subset(gene.sum, cnv)
    obs.epi = sum(exon.cnv.genes %in% epilepsy.genes)
    exp.epi.size = mclapply(1:nb.perm, function(ii) {
        gene.cont.ii = lapply(levels(gene.sum.cnv$size.class), function(sc) {
            sample(which(gene.sum$size.class == sc), sum(gene.sum.cnv$size.class == 
                sc))
        })
        gene.cont = gene.sum$geneName[unlist(gene.cont.ii)]
        sum(gene.cont %in% epilepsy.genes)
    }, mc.cores = NB.CORES)
    list(exp = unlist(exp.epi.size), obs = obs.epi, sum.df = data.frame(fold.enr = obs.epi/mean(unlist(exp.epi.size)), 
        gene = nrow(gene.sum.cnv), gene.epi = obs.epi, gene.epi.cont = mean(unlist(exp.epi.size)), 
        pv = (sum(obs.epi <= unlist(exp.epi.size)) + 1)/(length(exp.epi.size) + 
            1)))
}
```

### Exonic deletions absent from the public databases

``` r
test.epi.genes = list(allGenesDelNoDB = testGenesCnvs(cnv.all %>% filter(project == 
    "patients", z < 0, prop.db == 0), nb.perm = 1000, nb.class = 1), sizeGenesDelNoDB = testGenesCnvs(cnv.all %>% 
    filter(project == "patients", z < 0, prop.db == 0), nb.perm = 1000, nb.class = 5))

test.epi.genes.sum = do.call(rbind, lapply(names(test.epi.genes), function(x) data.frame(test = x, 
    test.epi.genes[[x]]$sum.df)))
kable(test.epi.genes.sum)
```

| test             |  fold.enr|  gene|  gene.epi|  gene.epi.cont|        pv|
|:-----------------|---------:|-----:|---------:|--------------:|---------:|
| allGenesDelNoDB  |  2.642624|   921|        17|          6.433|  0.000999|
| sizeGenesDelNoDB |  1.764217|   921|        17|          9.636|  0.018981|

``` r
test.df = rbind(data.frame(exp = test.epi.genes$allGenesDelNoDB$exp, test = "all genes"), 
    data.frame(exp = test.epi.genes$sizeGenesDelNoDB$exp, test = "size-controlled genes"))
ggplot(test.df, aes(x = exp, fill = test)) + geom_histogram(position = "dodge", 
    binwidth = 1) + theme_bw() + geom_vline(xintercept = test.epi.genes$sizeGenesDelNoDB$obs, 
    linetype = 2) + xlab("number of epilepsy genes among sampled genes") + ylab("number of sampling") + 
    scale_fill_brewer(palette = "Set2", name = "sampling") + theme(legend.position = "bottom") + 
    ggtitle("Genes hit deletions never seen in public databases") + annotate("text", 
    x = test.epi.genes$sizeGenesDelNoDB$obs, y = 100, label = "observed", vjust = -1, 
    angle = -90)
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-16-1.png)

### Exploring the enrichments

Here again, the number of permutations has been decreased compared to what was used in the study.

``` r
cnv.sets = cnv.all %>% mutate(type = ifelse(z < 0, "deletion", "duplication")) %>% 
    select(project, sample, chr, start, end, type) %>% as.data.frame
cnv.sets$prop.db = cnv.sets %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    dbProp(., svs.gr)

load("../data/cnvs-PopSV-twin-5kbp-FDR001.RData")
tw.sets = res.df %>% mutate(project = "Twins", type = ifelse(z < 0, "deletion", 
    "duplication")) %>% select(project, sample, chr, start, end, type)
tw.sets$prop.db = tw.sets %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>% 
    dbProp(., subset(svs.gr, project != "PopSV"))
cnv.sets = rbind(cnv.sets, tw.sets)
```

``` r
freq.r = c(1, 0.1, 0.01, 0.001, 1e-04, 0)
enr.epi.freq = lapply(freq.r, function(freq.max) {
    cnv.sets %>% filter(prop.db <= freq.max) %>% group_by(project, type) %>% 
        do(testGenesCnvs(., nb.perm = 1000)$sum.df) %>% mutate(freq = freq.max)
})
enr.epi.freq = do.call(rbind, enr.epi.freq)
```

``` r
ggplot(enr.epi.freq, aes(x = factor(freq), y = fold.enr, colour = project, group = project)) + 
    geom_line() + geom_point(aes(size = cut(pv, c(0, 0.05, 1)))) + theme_bw() + 
    ylim(0, max(enr.epi.freq$fold.enr)) + geom_hline(yintercept = 1, linetype = 3) + 
    scale_size_manual(name = "P-value", labels = c("<0.05", ">0.05"), values = c(4, 
        2)) + xlab("CNV frequency in public databases") + facet_grid(type ~ 
    .) + ylab("fold-enrichment") + scale_colour_brewer(name = "", palette = "Set1")
```

![](epilepsy-enrichmentPatterns_files/figure-markdown_github/unnamed-chunk-19-1.png)
