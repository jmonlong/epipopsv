PopSV - Methods benchmark using the CageKid dataset
===================================================

The [CageKid consortium](https://www.cng.fr/cagekid/) provides WGS for normal/tumor pairs of clear-cell renal carcinoma. Here we use the germline calls to evaluate the performance of the different CNV detection methods. For the vast majority of CNVs, we expect a germline variant to be present in the tumor.

Load packages, functions and data
---------------------------------

``` r
library(dplyr)
library(magrittr)
library(ggplot2)
library(PopSV)
library(GenomicRanges)
source("EpiPopSV-scripts.R")

## Get normal-tumor information
load("../data/cagekid-5kbp-files.RData")
files.df = files.df[sample.int(nrow(files.df), 45), ]
normals = subset(files.df, status == "normal")$sample
inds = subset(files.df, sample %in% normals)$individual
tumors = subset(files.df, individual %in% inds & status == "tumor")$sample

## CNVs from PopSV, FREEC, CNVnator, cn.MOPS, LUMPY
load("../data/cnvs-PopSV-cagekid-5kbp-FDR001.RData")
res.df$sig = res.df$qv
cnv.s = subset(res.df, sample %in% normals)
rm(res.df)
cnv.s = data.frame(method = "PopSV", set = "stringent", cnv.s[, c("sample", 
    "chr", "start", "end", "sig")])
load("../data/cnvs-PopSV-cagekid-5kbp-FDR05.RData")
res.df$sig = res.df$qv
cnv.l = subset(res.df, sample %in% tumors)
rm(res.df)
cnv.l = data.frame(method = "PopSV", set = "loose", cnv.l[, c("sample", "chr", 
    "start", "end", "sig")])
load("../data/cnvs-otherMethods-cagekid-5kbp.RData")
cnv.s = rbind(cnv.s, subset(others.df, sample %in% normals & set == "stringent")[, 
    c("method", "set", "sample", "chr", "start", "end", "sig")])
cnv.l = rbind(cnv.l, subset(others.df, sample %in% tumors & set == "loose")[, 
    c("method", "set", "sample", "chr", "start", "end", "sig")])
rm(others.df)

## Palette and method order
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
    "#CC79A7")
methods.f = c("LUMPY", "CNVnator", "cn.MOPS", "FREEC", "PopSV")
cnv.s$method = factor(as.character(cnv.s$method), levels = methods.f)
cnv.l$method = factor(as.character(cnv.l$method), levels = methods.f)
```

Systematic calls ?
------------------

How many systematic calls do we get in a typical sample ?

The distribution shows the average proportion of calls in one sample, grouped by their frequency in the full cohort.

``` r
cnv.s = cnv.s %>% group_by(method) %>% do(freq.range(., annotate.only = TRUE))
cnv.ave.samp = cnv.s %>% mutate(prop = cut(prop, seq(0, 1, 0.05), labels = seq(0.05, 
    1, 0.05))) %>% group_by(method, sample, prop) %>% summarize(call = n()) %>% 
    group_by(method, sample) %>% mutate(call = call/sum(call))
cnv.ave.samp %>% group_by(method, prop) %>% summarize(call = mean(call)) %>% 
    ggplot(aes(x = prop, y = call, fill = method)) + geom_bar(stat = "identity") + 
    facet_grid(method ~ ., scales = "free") + theme_bw() + scale_fill_manual(values = cbPalette) + 
    guides(fill = FALSE) + xlab("frequency") + ylab("proportion of calls")
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
cnv.ave.samp %>% filter(prop == "1") %>% ggplot(aes(x = method, y = call, fill = method)) + 
    geom_boxplot() + theme_bw() + coord_flip() + scale_fill_manual(values = cbPalette) + 
    guides(fill = FALSE) + ylim(0, 1) + ylab("proportion of systematic calls (>95% of the cohort)")
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-2-2.png)

Replication in the paired tumor
-------------------------------

``` r
samp.info = files.df[, c("sample", "individual")]
cnv.s %<>% filter(prop < 0.5)
cnv.s = merge(cnv.s, samp.info)
cnv.l = merge(cnv.l, samp.info)

concordance.nt <- function(cnv.df, cnv.2.df) {
    cnv.df$conc = overlapsAny(makeGRangesFromDataFrame(cnv.df), makeGRangesFromDataFrame(cnv.2.df))
    cnv.df
}

cnv.s = cnv.s %>% group_by(method) %>% do({
    subset(., individual %in% unique(subset(cnv.l, method == .$method[1])$individual))
})
cnv.s = cnv.s %>% group_by(method, individual) %>% do(concordance.nt(., subset(cnv.l, 
    method == .$method[1] & individual == .$individual[1]))) %>% ungroup
```

``` r
conc.nt = cnv.s %>% group_by(sample, method) %>% summarize(nb.c = sum(conc), 
    prop.c = mean(conc)) %>% mutate(method = factor(as.character(method), levels = methods.f)) %>% 
    group_by(method) %>% mutate(nb.c = winsorF(nb.c, med.u = 2))
ggplot(conc.nt, aes(x = method, y = prop.c)) + geom_boxplot(aes(fill = method)) + 
    theme_bw() + xlab("") + ylab("proportion of replicated calls per sample") + 
    ylim(0, 1) + coord_flip() + guides(fill = FALSE) + scale_fill_manual(values = cbPalette)
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
ggplot(conc.nt, aes(x = method, y = nb.c)) + geom_boxplot(aes(fill = method)) + 
    theme_bw() + xlab("") + ylab("number of replicated calls per sample") + 
    coord_flip() + guides(fill = FALSE) + scale_fill_manual(values = cbPalette)
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-4-2.png)

For LUMPY, CNVnator and PopSV we can further play with post-calling metrics. How does the performance change when playing with the calls' significance. For LUMPY, we use the number of supporting reads, for CNVnator the minimum of `eval1` and `eval2` and for PopSV the Q-value. More information on the `formatCalls-CNVnator-FREEC-LUMPY-cnMOPS.R` script (in the `others` folder).

``` r
conc.sum.sig <- function(df) {
    res = lapply(seq(0, 1, 0.05), function(sigq) {
        sigq.v = quantile(df$sig, probs = sigq)
        df %>% filter(sig < sigq.v) %>% group_by(sample) %>% summarize(nb.c = sum(conc), 
            prop.c = mean(conc)) %>% mutate(sigq = sigq)
    })
    do.call(rbind, res)
}
conc.tw.sig = cnv.s %>% filter(!is.na(sig)) %>% group_by(method) %>% do(conc.sum.sig(.))

conc.tw.sig %>% group_by(method, sigq) %>% summarize(nb.c = mean(nb.c), prop.u = quantile(prop.c, 
    1), prop.l = quantile(prop.c, 0), prop.c = mean(prop.c)) %>% ggplot(aes(x = nb.c, 
    y = prop.c, colour = method)) + geom_point() + theme_bw() + geom_line() + 
    geom_errorbar(aes(ymin = prop.l, ymax = prop.u)) + ylab("proportion of replicated calls per sample") + 
    xlab("number of replicated calls per sample") + scale_colour_manual(name = "", 
    values = cbPalette[c(1, 2, 5)]) + theme(legend.position = c(0.99, 0.01), 
    legend.justification = c(1, 0))
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-5-1.png)

Variants called by several methods vs unique to one method
----------------------------------------------------------

How many of the calls are unique to PopSV ? What is the proportion of calls from other methods that were found by PopSV ?

``` r
olMethods <- function(df) {
    ol = findOverlaps(makeGRangesFromDataFrame(df), makeGRangesFromDataFrame(df))
    ol.meth = tapply(df$method[subjectHits(ol)], queryHits(ol), function(x) unique(x))
    df = df[rep(as.numeric(names(ol.meth)), unlist(lapply(ol.meth, length))), 
        ]
    df$method2 = as.character(unlist(ol.meth))
    df
}
meth.df = cnv.s %>% group_by(sample) %>% do(olMethods(.))

call.samp = meth.df %>% group_by(method, sample, chr, start, end) %>% summarize(nb.meth = n()) %>% 
    group_by(sample, method) %>% summarize(tot.call = n())
meth.2d = meth.df %>% group_by(sample, method, method2) %>% summarize(nb.call = n()) %>% 
    merge(call.samp) %>% mutate(prop.call = nb.call/tot.call) %>% ungroup %>% 
    mutate(method = factor(method, levels = rev(methods.f)), method2 = factor(method2, 
        levels = methods.f))

meth.2d %>% group_by(method, method2) %>% summarize(prop.call = median(prop.call)) %>% 
    ggplot(aes(x = method, y = method2, fill = prop.call)) + geom_bin2d() + 
    theme_bw() + scale_fill_gradientn(name = "proportion\nof calls", colors = rev(terrain.colors(10))) + 
    xlab("calls from") + ylab("found by")
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
ggplot(subset(meth.2d, method != method2), aes(x = method2, fill = method, y = prop.call)) + 
    geom_bar(stat = "identity", position = "dodge", data = meth.2d %>% filter(method != 
        method2) %>% group_by(method, method2) %>% summarize(prop.call = median(prop.call))) + 
    geom_boxplot(alpha = 0, position = position_dodge(0.9)) + theme_bw() + scale_fill_manual(values = cbPalette, 
    name = "call from") + xlab("found by") + ylab("proportion of calls") + ylim(0, 
    1)
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-6-2.png)

If we focus only on regions called in at least two methods (out of 5).

``` r
call.samp.p = meth.df %>% group_by(method, sample, chr, start, end) %>% summarize(nb.meth = n()) %>% 
    filter(nb.meth > 1) %>% group_by(sample, method) %>% summarize(tot.call = n())
meth.2d.p = meth.df %>% group_by(method, sample, chr, start, end) %>% mutate(nb.meth = n()) %>% 
    filter(nb.meth > 1) %>% group_by(sample, method, method2) %>% summarize(nb.call = n()) %>% 
    merge(call.samp.p) %>% mutate(prop.call = nb.call/tot.call) %>% ungroup %>% 
    mutate(method = factor(method, levels = rev(methods.f)), method2 = factor(method2, 
        levels = methods.f))

meth.2d.p %>% group_by(method, method2) %>% summarize(prop.call = median(prop.call)) %>% 
    ggplot(aes(x = method, y = method2, fill = prop.call)) + geom_bin2d() + 
    theme_bw() + scale_fill_gradientn(name = "proportion\nof calls", colors = rev(terrain.colors(10))) + 
    xlab("calls from") + ylab("found by") + ggtitle("Calls in 2 methods or more")
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
ggplot(subset(meth.2d.p, method != method2), aes(x = method2, fill = method, 
    y = prop.call)) + geom_bar(stat = "identity", position = "dodge", data = meth.2d.p %>% 
    filter(method != method2) %>% group_by(method, method2) %>% summarize(prop.call = median(prop.call))) + 
    geom_boxplot(alpha = 0, position = position_dodge(0.9)) + theme_bw() + scale_fill_manual(values = cbPalette, 
    name = "call from") + xlab("found by") + ylab("proportion of calls") + ylim(0, 
    1) + ggtitle("Calls in 2 methods or more")
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-7-2.png)

Bias ?
------

Checking for potential bias when fragmented calls in one method versus stitched calls in another might duplicate the number of "calls" in the first method.

Let's see how many calls in one method overlaps several calls in another method

``` r
event.duplication.check <- function(cnv.o, methods = c("PopSV", "FREEC")) {
    gr1 = with(subset(cnv.o, method == methods[1]), GRanges(chr, IRanges(start, 
        end)))
    gr2 = with(subset(cnv.o, method == methods[2]), GRanges(chr, IRanges(start, 
        end)))
    ol = suppressWarnings(findOverlaps(gr1, gr2))
    if (length(ol) == 0) 
        return(data.frame(nb.ol = NA, count = NA, method = NA))
    ol.1 = as.data.frame(table(table(queryHits(ol))))
    ol.1$method = methods[1]
    ol.2 = as.data.frame(table(table(subjectHits(ol))))
    ol.2$method = methods[2]
    ol.df = rbind(ol.1, ol.2)
    colnames(ol.df)[1:2] = c("nb.ol", "count")
    ol.df$nb.ol = as.integer(as.character(ol.df$nb.ol))
    ol.df
}
dup.check = cnv.s %>% group_by(sample) %>% do({
    other.meth = setdiff(unique(.$method), "PopSV")
    tobind = lapply(other.meth, function(meth) {
        data.frame(comp = paste0("PopSV-", meth), event.duplication.check(., 
            methods = c("PopSV", meth)))
    })
    do.call(rbind, tobind)
})
dup.check = as.data.frame(dup.check)
dup.check$method = factor(as.character(dup.check$method), levels = methods.f)
```

``` r
ggplot(subset(dup.check, !is.na(nb.ol)), aes(x = winsorF(nb.ol, 3), y = count, 
    fill = method, group = paste(method, winsorF(nb.ol, 3)))) + geom_boxplot() + 
    theme_bw() + xlab("overlapping X calls from other method") + ylab("number of calls per sample") + 
    facet_wrap(~comp, scales = "free") + theme(text = element_text(size = 18)) + 
    scale_fill_manual(values = cbPalette) + scale_x_continuous(breaks = 1:3, 
    labels = c(1, 2, "3+"))
```

![](PopSV-methodBenchmark-cagekid_files/figure-markdown_github/unnamed-chunk-9-1.png)
