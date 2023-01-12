
library(vioplot)

if (!('d' %in% ls())){
    d <- read.table('counts_fig2B.txt', header=TRUE)
}

nb.simu <- 1000

sample.sizes <- c(22,21,24,19)
names(sample.sizes) <- c('CTA','LPO','TUR','TJE')
maf <- 0.1

pdf(sprintf('fig2B_MAF_%.2f.pdf', maf),  10, 7, title='Expected homozygote LOF variants')
par(mfrow=c(1,4))


pops <- unique(d$pop)
for (pop in pops){
    cat(pop,'\n')
    tmp.d <- d[d$pop==pop,]
    max.nb.alt <- maf * (2*sample.sizes[pop])
    tmp.d <- tmp.d[tmp.d$alt <= max.nb.alt,] 
    
    syns <- tmp.d[tmp.d$annot=='VEP_SS',]
    lof <- tmp.d[tmp.d$annot=='VEP_LOF',]
    expected <- c()
    ## we extract the frequencies (# alt allele) of LOF variants
    counts.nb.alt <- table(lof$alt)
    freqs <- as.numeric(names(counts.nb.alt))
    ## for each simu
    for (i in 1:nb.simu){
        random <- c()
        for (j in 1:length(counts.nb.alt)){ # for each frequency class
            tmp.freq <- freqs[j]
            nb.wanted <- counts.nb.alt[j]
            ## we identify the SYN variants with same frequency
            to.sample.from <- syns[syns$alt==tmp.freq,]
            ## we randomly draw from them
            tmp <- sample(to.sample.from$hom, nb.wanted)
            ## aggregate result
            random <- c(random, length(which(tmp > 0)))
        }
        expected <- c(expected, sum(random))
    }
    observed <- length(which(lof$hom > 0))
    vioplot(expected, ylim=range(c(expected, observed)), main=sprintf('%s MAF %.1f', pop, maf), xlab='', ylab='# variants with homozygote genotypes',xaxt='n')
    quantiles <- quantile(expected, probs=c(.025,.975))
    abline(h=quantiles,lty=2)
    points(1, observed, col='red',pch=19)
}

## max.freq <- .1
## for (pop in pops){
##     cat(pop,'\n')
##     tmp <- system(sprintf('wc -l ../samples_%s.txt', pop), intern=TRUE)
##     nb.sample <- as.integer(strsplit(tmp,' ')[[1]][1])
##     tmp.d <- d[d$pop==pop,]
##     syns <- tmp.d[tmp.d$annot=='VEP_SS',]
##     lof <- tmp.d[tmp.d$annot=='VEP_LOF',]
##     expected <- c()
##     ## we extract the frequencies (# alt allele) of LOF variants
##     counts.nb.alt <- table(lof$alt)
##     freqs <- as.numeric(names(counts.nb.alt))
##     ## only keep freqs less than what we want
##     to.keep <- freqs / (2*nb.sample) < max.freq
##     counts.nb.alt <- counts.nb.alt[to.keep]
##     freqs <- as.numeric(names(counts.nb.alt))
##     ## for each simu
##     for (i in 1:nb.simu){
##         random <- c()
##         for (j in 1:length(counts.nb.alt)){ # for each frequency class
##             tmp.freq <- freqs[j]
##             nb.wanted <- counts.nb.alt[j]
##             ## we identify the SYN variants with same frequency
##             to.sample.from <- syns[syns$alt==tmp.freq,]
##             ## we randomly draw from them
##             tmp <- sample(to.sample.from$hom, nb.wanted)
##             ## aggregate result
##             random <- c(random, length(which(tmp > 0)))
##         }
##         expected <- c(expected, sum(random))
##     }
##     ## compute statistic we are interested in
##     freqs <- lof$alt / (2*nb.sample)
##     to.keep <- freqs < max.freq
##     observed <- length(which(lof$hom[to.keep] > 0))

##     ## plot result
##     vioplot(expected, ylim=range(c(expected, observed)), main=pop, xlab='', ylab='# variants with homozygote genotypes',xaxt='n')
##     quantiles <- quantile(expected, probs=c(.025,.975))
##     abline(h=quantiles,lty=2)
##     points(1, observed, col='red',pch=19)
## }


dev.off()
