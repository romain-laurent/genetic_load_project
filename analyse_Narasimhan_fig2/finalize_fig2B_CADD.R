
library(vioplot)

if (!('d' %in% ls())){
    d <- read.table('counts_fig2B_CADD.txt', header=TRUE)
}

nb.simu <- 1000

sample.sizes <- c(22,21,24,19)
names(sample.sizes) <- c('CTA','LPO','TUR','TJE')
maf <- 1

pdf('fig2B_CADD.pdf',  10, 7, title='Expected homozygote LOF variants')
par(mfrow=c(1,4))


pops <- unique(d$pop)
for (pop in pops){
    cat(pop,'\n')
    tmp.d <- d[d$pop==pop,]
    if (maf < 1){
        max.nb.alt <- maf * (2*sample.sizes[pop])
        tmp.d <- tmp.d[tmp.d$alt <= max.nb.alt,]
    }
    
    syns <- tmp.d[tmp.d$annot=='CADD_inf_10',]
    lof <- tmp.d[tmp.d$annot=='CADD_sup_28',]
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


dev.off()
