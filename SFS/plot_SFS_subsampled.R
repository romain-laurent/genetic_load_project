
pops <- c('CTA','LPO','TUR','TJE')
sample.size <- 19
possible.counts <- 1:(2*sample.size-1)

## read all data once
if (!('counts' %in% ls())){
    
    ## prepare datastruct for results
    counts <- matrix(NA, nrow=length(pops), ncol=2*sample.size-1)
    rownames(counts) <- pops
    colnames(counts) <- possible.counts

    ## read allele counts
    for (pop in pops){
        if (!(file.exists(sprintf('%s_SFS_subsampled.Rdata',pop)))){
            cat(pop,'\n')
            tmp <- read.table(sprintf('SFS_%s_subsampled.frq.count', pop), skip=1)
            save(tmp, file=sprintf('%s_SFS_subsampled.Rdata',pop))
        } else {
            cat('Reading',pop,'from Rdata...\n')
            load(sprintf('%s_SFS.Rdata',pop))
        }
        ## remove monomorphic markers, i.e. derived counts == 0 or reference counts == 0
        colnames(tmp) <- c('chrom','pos','num.alleles','num.chrom','ref','alt')
        tmp <- tmp[tmp$alt != 0,]
        tmp <- tmp[tmp$ref != 0,]
        ## make alt column a factor with levels = possible counts (useful to have a table with possibly missing categories)
        tmp$alt <- factor(tmp$alt, levels=possible.counts)
        ## compute SFS + save result
        counts[pop,] <- table(tmp$alt)/nrow(tmp)
    }
}

## plot result
pdf('SFS_subsampled.pdf', 15, 10, title='SFS')
cols <- c('tomato1', 'gold', 'seagreen3', 'lightskyblue1')
names(cols) <- pops
barplot(counts, beside=TRUE, col=cols, xlab='Derived allele count', ylab='Fraction of sites', legend.text=pops)
dev.off()
