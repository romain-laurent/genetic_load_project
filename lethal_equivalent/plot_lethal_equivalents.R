
plot.bootstrap <- function(neutral, non.neutral){
    filename <- sprintf('bootstrap_%s_%s.txt', neutral, non.neutral)
    d <- read.table(filename, header=TRUE)
    
    real.values <- d[d$idx==0,2:5]
    boot.values <- d[d$idx > 0,2:5]
    
    ic <- apply(boot.values, 2, quantile, probs=c(0.025,.975))
    
    cols <- c('tomato1', 'gold', 'seagreen3', 'lightskyblue1')
    ylims <- range(ic)
    plot.main <- sprintf('%s vs %s', neutral, non.neutral)
    plot(1:4, as.numeric(real.values), col=cols, pch=19, xaxt='n', xlab='', ylab='# variants', ylim=ylims, xlim=c(.75,4.25), main=plot.main)
    arrows(x0=1:4, y0=ic[1,], y1=ic[2,], col=cols, code=3, angle=90, lwd=2)
    axis(1, at=1:4, labels=colnames(boot.values))
    abline(h=0, col='red', lwd=2, lty=2)
}

pdf('bootstrap_lethal_equivalents.pdf', 15, 15, title='Boostrap lethal equivalents')
par(mfrow=c(2,2))
plot.bootstrap('VEP_SS','VEP_LOF')
plot.bootstrap('CADD_inf_10', 'CADD_sup_28')
plot.bootstrap('CADD_inf_10', 'CADD_15_28')
plot.bootstrap('CADD_inf_10', 'CADD_10_15')

dev.off()
