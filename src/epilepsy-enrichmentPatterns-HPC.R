library(BatchJobs)

permExCNV.f <- function(cnv.prof){
  library(PopSV)
  library(dplyr)
  library(magrittr)
  library(GenomicRanges)
  library(parallel)
  library(tidyr)
  ## Load functions and genomic annotations
  source("EpiPopSV-scripts.R")
  load("genomicFeatures.RData")
  cnv.all = read.table("cnvs-PopSV-Epilepsy-198affected-301controls-5kb.tsv.gz", header=TRUE, as.is=TRUE, sep="\t")
  exac = read.csv("fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", as.is=TRUE, sep="\t")
  ## Exons
  exons.grl = list(exon=subset(gen.feat.l$exon, geneType=="protein_coding"))
  ## ExAC lof intolerant
  gene.intolerant = subset(exac, pLI>.9)$gene
  exons.grl$exon.pli = subset(gen.feat.l$exon, geneName %in% gene.intolerant)
  ## Convert to GRanges
  exons.grl = lapply(exons.grl, makeGRangesFromDataFrame, keep.extra.columns=TRUE)
  ## Load and annotate CNVs
  cnv.all = cnv.all %>% mutate(project=ifelse(project=="affected", "patients", "controls"))
  cnv.all.gr = makeGRangesFromDataFrame(cnv.all, keep.extra.columns=TRUE)
  cnv.all$exon = overlapsAny(cnv.all.gr, exons.grl$exon)
  info.df = cnv.all %>% select(sample, project) %>% unique
  sampToProj = info.df$project
  names(sampToProj) = info.df$sample
  samps = info.df$sample
  permF <-  function(rr, cnv.df, feat.grl, nperms=1000, ss.nb=150){
    message(rr)
    info.df %<>% group_by(project) %>% filter(sample %in% sample(unique(sample), ss.nb))
    cat.all = cnv.df %>% filter(sample %in% info.df$sample) %>% group_by(project) %>% do(reduceDf(.)) %>% ungroup %>% mutate(sample=1:n(), control=FALSE) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    idsToProj = cat.all$project
    cat.null = draw.controls(cat.all, list(centel=gen.feat.l$centel), nb.cores=1)
    cat.null$project = idsToProj[cat.null$sample]
    cat.null$control = TRUE
    mcols(cat.null) = mcols(cat.null)[,c('project','sample','control')]
    cat.all = c(cat.all, cat.null)
    res = lapply(1:length(feat.grl), function(ii){
      cat.all$ol = overlapsAny(cat.all, feat.grl[[ii]])
      obs = mcols(cat.all) %>% as.data.frame %>% group_by(control, project) %>% summarize(nb=sum(ol), prop=(1+nb)/(1+n())) %>% select(-nb) %>% group_by(project) %>% spread(control, prop, sep='.') %>% mutate(feat=names(feat.grl)[ii], rep=rr, set="obs", perm=NA)
      exp = lapply(1:nperms, function(perm){
        idsToProj = sample(idsToProj)
        mcols(cat.all) %>% as.data.frame %>% mutate(project=idsToProj[sample]) %>% group_by(control, project) %>% summarize(nb=sum(ol), prop=(1+nb)/(1+n())) %>% select(-nb) %>% group_by(project) %>% spread(control, prop, sep='.') %>% mutate(feat=names(feat.grl)[ii], rep=rr, set="exp", perm=perm)
      })
      exp = do.call(rbind, exp)
      rbind(obs, exp)
    })
    do.call(rbind, res)
  }
  NB.CORES = 10
  NB.SUBSAMP = 100
  NB.PERMS = 10000
  cnv.size = gsub('(.*)-.*-.*', '\\1', cnv.prof)
  cnv.freq = gsub('.*-(.*)', '\\1', cnv.prof)
  if(cnv.freq == 'rare'){
    load("SVdatabase.RData")
    if(cnv.size=='small'){
      cnv.all$prop.db = dbProp(cnv.all.gr, svs.gr)
      cnv.all = filter(cnv.all, prop.db<.01)
    } else {
      cnv.all$prop.db = dbProp(cnv.all.gr, svs.gr, min.db.span = 0.5) # 50% overlap
      cnv.all = filter(cnv.all, prop.db<.01)
    }
  } 
  if(cnv.size == 'small'){
    cnv.all = filter(cnv.all, (end-start+1) < 5e4)
  }
  if(cnv.size == 'large'){
    cnv.all = filter(cnv.all, (end-start+1) >= 5e4)
  }
  res = do.call(rbind, mclapply(1:NB.SUBSAMP, permF, cnv.df=cnv.all, feat.grl=exons.grl, nperms=NB.PERMS, mc.cores=NB.CORES))
  enr.obs = res %>% filter(is.na(perm)) %>% group_by(feat,project, rep) %>% summarize(enr=control.FALSE/control.TRUE) %>% mutate(freq=cnv.freq, size=cnv.size)
  enr.exp.diff = res %>% filter(!is.na(perm)) %>% group_by(perm,feat,project,rep) %>% summarize(enr=control.FALSE/control.TRUE) %>% group_by(perm,feat,project) %>% summarize(enr=median(enr)) %>% group_by(perm, feat) %>% summarize(diff=enr[project=='patients']-enr[project=='controls']) %>% mutate(freq=cnv.freq, size=cnv.size)
  list(enr.obs=enr.obs, enr.exp.diff=enr.exp.diff)
}

cnv.profs = c('small-all', 'large-all', 'small-rare', 'large-rare')
permExCNV.reg <- makeRegistry(id="permExCNV")
batchMap(permExCNV.reg, permExCNV.f, cnv.profs)
submitJobs(permExCNV.reg, findNotSubmitted(permExCNV.reg), resources=list(walltime="12:0:0", nodes="1", cores="12", supervisor.group="bws-221-ae", queue='metaq'))
showStatus(permExCNV.reg)

permExCnv = reduceResultsList(permExCNV.reg)
save(permExCnv, file='../data/permExCnv.RData')

