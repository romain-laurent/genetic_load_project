
pops <- c('CTA','LPO','TUR','TJE')

## find plot limits
bp.max <- 0
cM.max <- 0
for (pop in pops){
    d.hbd <- read.table(sprintf('results/flod-%s/HBD.segments.txt', pop), header=TRUE)
    d.hbd$len.bp <- d.hbd$FIN_bp - d.hbd$DEB_bp
    d.hbd$len.cM <- d.hbd$FIN_cM - d.hbd$DEB_cM
    for (indiv in unique(d.hbd$IID)){
        tmp.d <- d.hbd[d.hbd$IID == indiv,]
        length.bp <- sum(tmp.d$len.bp)
        length.cM <- sum(tmp.d$len.cM)
        length.bp  <- length.bp / 1e6
        bp.max <- max(bp.max, length.bp)
        cM.max <- max(cM.max, length.cM)
    }
    
}

all.indivs.hbd <- NULL
pdf('FEstim_HBD_segment_per_indiv.pdf', 15, 7, title='HBD FEstim per indiv')
par(mfrow=c(2,2))
for (pop in pops){
    
    d.festim <- read.table(sprintf('results/festim-%s.summary', pop), header=TRUE)
    
    ## extract inferred mating types for parents
    mat.type <- d.festim[,14:ncol(d.festim)]
    mat.type <- colnames(mat.type)[apply(mat.type, 1, which.max)]
    mat.type <- sub('X', '', mat.type)
    names(mat.type) <- d.festim$IID
    Fs <- d.festim$F_MEDIAN
    names(Fs) <- d.festim$IID
    
    ## now we work on HBD segments
    d.hbd <- read.table(sprintf('results/flod-%s/HBD.segments.txt', pop), header=TRUE)
    d.hbd$len.bp <- d.hbd$FIN_bp - d.hbd$DEB_bp
    d.hbd$len.cM <- d.hbd$FIN_cM - d.hbd$DEB_cM
    ## compute total HBD length per individual
    final.d.hbd <- NULL
    for (indiv in names(mat.type)){
        tmp.d <- d.hbd[d.hbd$IID == indiv,]
        
        length.bp <- 0
        length.cM <- 0
        if (nrow(tmp.d) > 0){
            length.bp <- sum(tmp.d$len.bp)
            length.cM <- sum(tmp.d$len.cM)
        }
        tmp <- data.frame(indiv=indiv, bp=length.bp, cM=length.cM)
        final.d.hbd <- rbind(final.d.hbd, tmp)
    }
    
    ## plot
    final.d.hbd$bp <- final.d.hbd$bp / 1e6
    cols <- cbind(c('1C','2C','2x1C','AV','OUT'),c('red','blue','yellow','brown','black'))
    my.cols <- cols[match(mat.type, cols[,1]), 2]
    xs <- barplot(bp~indiv, data=final.d.hbd, las=2, col=my.cols, xlab='', main=pop, ylim=c(0, bp.max), ylab='Total HBD segments length (Mb)')

    text(xs, final.d.hbd$bp, Fs, xpd=TRUE, srt=90, adj=c(-0.1,0.4))
    
    ## ## plot
    ## cols <- cbind(c('1C','2C','2x1C','AV','OUT'),c('red','blue','yellow','brown','black'))
    ## my.cols <- cols[match(mat.type, cols[,1]), 2]
    ## x <- barplot(cM~indiv, data=final.d.hbd, las=2, col=my.cols, xlab='', main=pop, ylim=c(0, cM.max))

    legend('topleft', legend=cols[,1], fill=cols[,2])


    all.indivs.hbd <- rbind(all.indivs.hbd, final.d.hbd)
}
dev.off()

pdf('FEstim_bp_cM_correlation.pdf', 10, 7, title='Correlation bp/cM')
tmp <- cor.test(all.indivs.hbd$bp, all.indivs.hbd$cM, method='spearman')
main <- sprintf('Spearman rho = %.3f, p = %.1e', tmp$estimate, tmp$p.value)
plot(all.indivs.hbd$bp, all.indivs.hbd$cM, xlab='Total HBD length (Mb)', ylab='Total HBD length (cM)', main=main)
abline(lm(cM~bp, data=all.indivs.hbd), col='red', lty=2)
dev.off()
