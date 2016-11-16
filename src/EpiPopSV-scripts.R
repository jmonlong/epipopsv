##
## Used for the epilepsy analysis
##

winsorF <- function(x, u=NULL){
    if(is.null(u)) u = min(median(x, na.rm=TRUE)*3, quantile(x, .95, na.rm=TRUE))
    if(any(x>u)) x[x>u] = u
    x
}
dbProp <- function(gr, db.gr, min.db.span=0){
    ol = findOverlaps(gr, db.gr)
    prop = rep(0,length(gr))
    ol.span = width(pintersect(gr[queryHits(ol)], db.gr[subjectHits(ol)]))/width(db.gr)[subjectHits(ol)]
    ol.span.ii = which(ol.span>min.db.span)
    ol.t = tapply(db.gr$prop[subjectHits(ol)[ol.span.ii]], queryHits(ol)[ol.span.ii], max)
    prop[as.numeric(names(ol.t))] = as.numeric(ol.t)
    prop
}
rmOL <- function(df, nb.ol=1){
    df.s = df %>% group_by(project, sample) %>% summarize(kb=sum((end-start)/1e3)) %>% arrange(desc(kb))
    print(head(as.data.frame(df.s), nb.ol))
    subset(df, !(sample %in% head(df.s$sample, nb.ol)))
}
decdf <- function(df, coln="exon.epi.d"){
    cnv.ecdf = ecdf(unlist(df[,coln]))
    res = data.frame(d=seq(0,3e5,200))
    res$cnb = cnv.ecdf(res$d) * nrow(df)
    res
}
freqSS <- function(df, sample.size=198, nb.rep=100, nb.cores=3){
    if(length(unique(df$sample)) > sample.size){
        res = mclapply(1:nb.rep, function(ii){
                           df %>% filter(sample %in% sample(unique(sample),sample.size)) %>% do(freq.range(.,annotate.only=TRUE)) %>% mutate(rep=ii) %>% as.data.frame
                       }, mc.cores=nb.cores)
        return(do.call(rbind, res))
    } else {
        return(df %>% freq.range(annotate.only=TRUE) %>% mutate(rep=1))
    }
}
reduceDf <- function(df, stitch.dist=0){
  df %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>% reduce(min.gapwidth=stitch.dist) %>% as.data.frame %>% mutate(chr=seqnames) %>% select(chr, start, end)
}
