
pops <- c('CTA','LPO','TUR','TJE')
sample.sizes <- c(22,21,24,19)
names(sample.sizes) <- pops

## read data only once
if (!('d' %in% ls())){
    d <- read.table('../analyse_Narasimhan_fig2/counts_fig2B.txt', header=TRUE)
}


neutral <- 'VEP_SS'
non.neutral <- 'VEP_LOF'
## for each pop
for (pop in pops){
#    cat(pop,'\n')
    ## only get data for this pop
    tmp.d <- d[d$pop==pop,]
    sample.size <- sample.sizes[pop]

    ## get possible alternate allele counts (i.e. allele frequencies)
    alt.counts <- unique(tmp.d$alt)
    ## for each allele frequency, count the depletion of homozygous
    depletion <- rep(NA, length(alt.counts))
    for (idx in 1:length(alt.counts)){
        alt.count <- alt.counts[idx]
        nb.hom.neutral <- mean(tmp.d$hom[tmp.d$annot==neutral & tmp.d$alt==alt.count])
        nb.hom.non.neutral <- mean(tmp.d$hom[tmp.d$annot==non.neutral & tmp.d$alt==alt.count])
        depletion[idx] <- nb.hom.neutral - nb.hom.non.neutral
    }
    depletion <- depletion * alt.counts / (2*sample.size)
    depletion <- depletion[!is.nan(depletion)]
    cat(pop, sum(depletion),'\n')
    
}
