SV databases
============

``` r
library(ggplot2)
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(data.table)
library(knitr)
library(tidyr)
library(PopSV)
```

Several studies looked at Structural Variation (SV) across many samples. Usually, they report the allele frequency. For simple deletions/insertions, the concept of allele and allele frequency is fine, but it can be more complex for other SVs, especially CNVs. In practice, we are more interested in the *proportion of individuals with a SV*. This is the information we would like to know when we observe a (recurrent) SV: how often an individual has a SV in this particular region ?

For this reason, we download the calls from public catalogs, compute the frequency as the *proportion of samples with a SV*, and combine them into a R object for future use.

1000 Genomes Project
--------------------

``` r
if (!file.exists("ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz")) download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz", 
    "ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz")
tgp = read.table("ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz", 
    as.is = TRUE, nrows = 10)
tgp.cc = c("character", "integer", rep("NULL", 5), "character", "NULL", rep("character", 
    ncol(tgp) - 9))
tgp = read.table("ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz", 
    as.is = TRUE, colClasses = tgp.cc)
colnames(tgp)[1:3] = c("chr", "start", "info")
samples = 4:ncol(tgp)
tgp$type = gsub(".*SVTYPE=([^;]+);*.*", "\\1", tgp$info)
tgp$end = as.numeric(gsub(".*END=([^;]+);*.*", "\\1", tgp$info))
tgp$end = ifelse(is.na(tgp$end), tgp$start, tgp$end)
tgp$prop = apply(tgp[, samples], 1, function(x) mean(x != "0|0" & x != "0"))
tgp$af = gsub(".*;AF=([^;]+);*.*", "\\1", tgp$info)
tgp$af = as.numeric(unlist(lapply(strsplit(tgp$af, ","), function(x) sum(as.numeric(x)))))
tgp$size = as.numeric(gsub(".*SVLEN=([^;]+);*.*", "\\1", tgp$info))
tgp$size = ifelse(is.na(tgp$size), tgp$end - tgp$start, tgp$size)
tgp %<>% mutate(project = "1KGP", type = ifelse(grepl("DEL", type), "DEL", type), 
    type = ifelse(type %in% c("ALU", "SVA", "LINE1"), "MEI", type)) %>% select(chr, 
    start, end, size, prop, af, type, project)
```

The calls are downloaded from [ALL.wgs.integrated\_sv\_map\_v2.20130502.svs.genotypes.vcf.gz](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz).

### Allele frequency vs sample proportion

The 1000 Genomes Project provides both the allele frequencies and the genotype for each sample. We can compare the usual allele frequency estimate with the one we want to use (proportion of samples). Assuming Hardy-Weinberg equilibrium, we expect (1-(1-AF)^2) of the samples to have a mutated allele.

``` r
hw.exp = data.frame(af = seq(0, 1, 0.01))
hw.exp$prop = 1 - (1 - hw.exp$af) * (1 - hw.exp$af)
ggplot(tgp, aes(x = af, y = prop)) + geom_point(aes(colour = chr %in% 1:22), 
    alpha = 0.5) + geom_line(data = hw.exp, linetype = 1) + geom_line(aes(y = 1.1 * 
    prop), data = hw.exp, linetype = 2) + facet_wrap(~type) + theme_bw() + scale_colour_hue(name = "autosome") + 
    xlab("allele frequency") + ylab("proportion of samples")
```

![](annotation-SVdatabases_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
tgp$prop.exp = 1 - (1 - tgp$af) * (1 - tgp$af)
```

Most of the variants follow the HB expectation (plain line). Many of the ones that deviates from HB affect **more samples than expected**. This is not so surprising, these variants were most likely depleted in homozygous variant by selective pressure. The few variants that affect less samples than expected are on chromosomes X and Y.

In practice, if we estimate the proportion of samples from the allele frequency we will **under-estimate** it. Here for example, it would be under-estimated for 83.4% of the variants. Instead we could use a slightly **skewed HB expectation** (dotted line). Then, the under-estimation would only affect 0.7% of the variants.

Conclusions:

-   If we can, it's more accurate to compute the *proportion of samples* directly from the calls.
-   If not, we could always use the allele frequency and the skewed HB formula, to be sure we over-estimate rather than under-estimate the proportion of samples.

``` r
tgp %<>% select(chr, start, end, size, prop, type, project)
```

GoNL
----

``` r
if (!file.exists("20160525_GoNL_AF_genotyped_SVs.vcf.gz")) download.file("https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release6/20160525_GoNL_AF_genotyped_SVs.vcf.gz", 
    "20160525_GoNL_AF_genotyped_SVs.vcf.gz")
gonl = read.table("20160525_GoNL_AF_genotyped_SVs.vcf.gz", as.is = TRUE, sep = "\t", 
    quote = "")
colnames(gonl) = c("chr", "start", "id", "ref", "alt", "qual", "filter", "info")
gonl$type = gsub(".*SVTYPE=([^;]+);*.*", "\\1", gonl$info)
gonl$af = as.numeric(gsub(".*;AF=([^;]+);*.*", "\\1", gonl$info))
gonl$size = abs(as.numeric(gsub(".*SVLEN=([^;]+);*.*", "\\1", gonl$info)))
gonl$end = as.numeric(gsub(".*END=([^;]+);*.*", "\\1", gonl$info))
gonl$end = ifelse(gonl$type == "DEL", gonl$start + gonl$size, gonl$end)
gonl$end = ifelse(is.na(gonl$end), gonl$start, gonl$end)
gonl$size = ifelse(is.na(gonl$size), gonl$end - gonl$start, gonl$size)
gonl %<>% filter(type != "TRA") %>% mutate(prop = 1.1 * (1 - (1 - af) * (1 - 
    af)), prop = ifesle(prop > 1, 1, prop), project = "GoNL", type = ifelse(type == 
    "INS:ME", "MEI", type)) %>% select(chr, start, end, size, prop, type, project)
```

GoNL calls are downloaded from [20160525\_GoNL\_AF\_genotyped\_SVs.vcf.gz](https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release6/20160525_GoNL_AF_genotyped_SVs.vcf.gz). They just provide the allele frequencies of the variants so we estimated the proportion of samples with an mutated allele. We used the conservative HB formula (see previously).

CNVs - Sudmant Science 2015
---------------------------

``` r
if (!file.exists("Sudmant-Science2015-S1.csv")) download.file("https://dl.dropboxusercontent.com/s/537ghamwr54tdxe/Sudmant-Science2015-S1.csv", 
    "Sudmant-Science2015-S1.csv")
sud = read.csv("Sudmant-Science2015-S1.csv", skip = 1)
samples = colnames(sud)[-(1:6)]
sud$nb = apply(sud[, samples], 1, function(x) sum(x != -1 & x != 2))
sud$prop = sud$nb/length(samples)
sud %<>% filter(nb > 0) %>% mutate(size = end - start, chr = gsub("chr", "", 
    contig), project = "SudmantScience2015") %>% select(chr, start, end, size, 
    prop, type, project)
```

From [Sudmant et al Science 2015](http://science.sciencemag.org/content/349/6253/aab3761), [Supplementary Table 1](http://science.sciencemag.org/highwire/filestream/633982/field_highwire_adjunct_files/1/aab3761_TableS1.xlsx) has all the CNV genotypes for all their samples.

Handsaker - Nature Genetics 2015
--------------------------------

``` r
if (!file.exists("1000G_phase1_cnv_genotypes_phased_25Jul2014.genotypes.vcf.gz")) download.file("http://www.broadinstitute.org/~handsake/mcnv_data/bulk/1000G_phase1_cnv_genotypes_phased_25Jul2014.genotypes.vcf.gz", 
    "1000G_phase1_cnv_genotypes_phased_25Jul2014.genotypes.vcf.gz", method = "wget")
hand = read.table("1000G_phase1_cnv_genotypes_phased_25Jul2014.genotypes.vcf.gz", 
    as.is = TRUE, sep = "\t", quote = "")
hand = hand[, c(1, 2, 8, 10:ncol(hand))]
hand = gather(hand, "sample", "gt", 4:ncol(hand))
colnames(hand)[1:3] = c("chr", "start", "info")
hand$type = gsub(".*SVTYPE=([^;]+);*.*", "\\1", hand$info)
hand$end = as.numeric(gsub(".*END=([^;]+);*.*", "\\1", hand$info))
hand$gt = unlist(lapply(strsplit(hand$gt, ":"), "[", 1))
hand = hand %<>% filter(gt != "0|0") %>% freq.range
hand %<>% mutate(type = "CNV", project = "GenomeSTRiP2", size = end - start) %>% 
    select(chr, start, end, size, prop, type, project)
```

From [Handsaker et al Nat. Genetics 2015](http://www.nature.com/ng/journal/v47/n3/full/ng.3200.html), GenomeSTRiP 2 calls are available in [1000G\_phase1\_cnv\_genotypes\_phased\_25Jul2014.genotypes.vcf.gz](http://www.broadinstitute.org/~handsake/mcnv_data/bulk/1000G_phase1_cnv_genotypes_phased_25Jul2014.genotypes.vcf.gz).

PopSV
-----

``` r
if (!file.exists("CNV-PopSV-Twin_CageKid_GoNL-germline.tsv")) download.file("https://ndownloader.figshare.com/files/3638574?private_link=ba79730bb87a1322480d", 
    "CNV-PopSV-Twin_CageKid_GoNL-germline.tsv")
popsv = read.table("CNV-PopSV-Twin_CageKid_GoNL-germline.tsv", header = TRUE, 
    as.is = TRUE)
popsv %<>% filter(!grepl("gonl", sample)) %>% select(chr, start, end, sample) %>% 
    freq.range %>% mutate(type = "CNV", project = "PopSV", size = end - start) %>% 
    select(chr, start, end, size, prop, type, project)
```

The catalog from PopSV is based on a Twin study and normals from CageKig project. The calls were fragmented and each sub-call is saved with the proportion of affected samples.

Merge
-----

All the CNVs are merged into one GRanges object, annotated with the frequency in each cohort.

``` r
svs.gr = do.call(rbind, list(popsv, tgp, gonl, sud, hand)) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
save(svs.gr, file = "../data/SVdatabase.RData")
```
