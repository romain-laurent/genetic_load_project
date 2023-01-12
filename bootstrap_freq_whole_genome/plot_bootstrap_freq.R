
comparisons <- c('CTA_LPO','TUR_TJE')

dtt <- read.table('result_bootstrap_freq_TUR_TJE.txt', header=TRUE)
dcl <- read.table('result_bootstrap_freq_CTA_LPO.txt', header=TRUE)

tmp <- rbind(dtt, dcl)
xlims <- range(tmp[,-1])

plot.conf.int <- function(values, d.values, probs=c(0.025, 0.975), col='black'){
    conf.int <- quantile(values, probs=probs)
    tmp.fun <- approxfun(d.values$x, d.values$y)

    ## first part
    final.x <- conf.int[1]
    final.y <- tmp.fun(final.x)
    xs <- d.values$x[d.values$x < final.x]
    ys <- d.values$y[d.values$x < final.x]
    xs <- c(xs, final.x, final.x)
    ys <- c(ys, final.y, 0)
    polygon(xs, ys, col=col, border=col)

    ## second part
    final.x <- conf.int[2]
    final.y <- tmp.fun(final.x)
    xs <- d.values$x[d.values$x > final.x]
    ys <- d.values$y[d.values$x > final.x]
    xs <- c(final.x, final.x, xs)
    ys <- c(0, final.y, ys)
    polygon(xs, ys, col=col, border=col)

}


pdf('bootstrap_freq_whole_genome.pdf', 15, 10, title='As much variants as possible')
par(mfrow=c(3,2))
for (idx.col in 2:5){
    real.dtt <- dtt[1,idx.col]
    real.dcl <- dcl[1,idx.col]
    tmp.dtt <- dtt[-1,idx.col]
    tmp.dcl <- dcl[-1,idx.col]

    d.dtt <- density(tmp.dtt)
    d.dcl <- density(tmp.dcl)
    ylims <- range(c(d.dtt$y, d.dcl$y))
    xlims <- range(c(d.dtt$x, d.dcl$x))
    
    plot(d.dtt, xlab='ratio', main=colnames(dtt)[idx.col], col='slategray3', lwd=2, xlim=xlims, ylim=ylims)
    lines(d.dcl, col='salmon1', lwd=1)
    abline(v=real.dtt, col='slategray3', lty=2, lwd=2)
    abline(v=real.dcl, col='salmon1', lty=2, lwd=1)

    ## plot 95% confidence intervals
    plot.conf.int(tmp.dtt, d.dtt, col='slategray3')
    plot.conf.int(tmp.dcl, d.dcl, col='salmon1')
    ## plot "theoretical value"
    abline(v=1, col='red')
    
    legend('topright', legend=c('CTA/LPO', 'TUR/TJE'), col=c('salmon1','slategray3'), lwd=c(1,2))
}

tmp.dtt <- dtt[-1,3] + dtt[-1,5]
tmp.dcl <- dcl[-1,3] + dcl[-1,5]
real.dtt <- dtt[1,3] + dtt[1,5]
real.dcl <- dcl[1,3] + dcl[1,5]


d.dtt <- density(tmp.dtt)
d.dcl <- density(tmp.dcl)
ylims <- range(c(d.dtt$y, d.dcl$y))
xlims <- range(c(d.dtt$x, d.dcl$x))

plot(d.dtt, xlab='ratio', main='sum all_der + all_anc', col='slategray3', lwd=2, xlim=xlims, ylim=ylims)
lines(d.dcl, col='salmon1', lwd=1)
abline(v=real.dtt, col='slategray3', lty=2, lwd=2)
abline(v=real.dcl, col='salmon1', lty=2, lwd=1)

## plot 95% confidence intervals
plot.conf.int(tmp.dtt, d.dtt, col='slategray3')
plot.conf.int(tmp.dcl, d.dcl, col='salmon1')
## plot "theoretical value"
abline(v=1, col='red')

legend('topright', legend=c('CTA/LPO', 'TUR/TJE'), col=c('salmon1','slategray3'), lwd=c(1,2))
dev.off()
