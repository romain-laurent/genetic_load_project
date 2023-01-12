
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


gerp.scores <- factor(c("GERP_inf_02", "GERP_02_04", "GERP_04_06", "GERP_sup_06"), levels=c("GERP_inf_02", "GERP_02_04", "GERP_04_06", "GERP_sup_06"), ordered=TRUE)
cadd.scores <- factor(c("CADD_inf_10", "CADD_10_15", "CADD_15_28", "CADD_sup_28"), levels=c("CADD_inf_10", "CADD_10_15", "CADD_15_28", "CADD_sup_28"), ordered=TRUE)
all.scores <- list(gerp.scores=gerp.scores, cadd.scores=cadd.scores)


dcl <- read.table('bootstrap_freq_CTA_LPO.txt', header=TRUE)
dtt <- read.table('bootstrap_freq_TUR_TJE.txt', header=TRUE)
dsc <- read.table('bootstrap_freq_SEA_CA.txt', header=TRUE)

pdf('bootstrap_freq_with_SEA_CA.pdf', 15, 10, title='Bootstrap freq')
for (i in 1:length(all.scores)){
    scores <- all.scores[[i]]
    par(mfrow=c(length(scores), 2))
    
    for (score.class in scores){
        tmp.dcl <- dcl[dcl$score.class==score.class,]
        tmp.dtt <- dtt[dtt$score.class==score.class,]
        tmp.dsc <- dsc[dsc$score.class==score.class,]

        # for all derived
        real.value.dcl <- tmp.dcl$all.der[tmp.dcl$idx == 0]
        real.value.dtt <- tmp.dtt$all.der[tmp.dtt$idx == 0]
        real.value.dsc <- tmp.dsc$all.der[tmp.dsc$idx == 0]
        boot.value.dcl <- tmp.dcl$all.der[tmp.dcl$idx != 0]
        boot.value.dtt <- tmp.dtt$all.der[tmp.dtt$idx != 0]
        boot.value.dsc <- tmp.dsc$all.der[tmp.dsc$idx != 0]

        density.dcl <- density(boot.value.dcl)
        density.dtt <- density(boot.value.dtt)
        density.dsc <- density(boot.value.dsc)
        xlims <- range(c(density.dcl$x, density.dtt$x, density.dsc$x))
        ylims <- range(c(density.dcl$y, density.dtt$y, density.dsc$y))

        plot(0,0,type='n', xlim=xlims, ylim=ylims, main=paste(score.class, 'All derived'), ylab='Density', xlab='Ratio')
        lines(density.dcl, col='salmon1')
        abline(v=real.value.dcl, col='salmon1', lty=2)
        lines(density.dtt, col='slategray3')
        abline(v=real.value.dtt, col='slategray3', lty=2)
        lines(density.dsc, col='green')
        abline(v=real.value.dsc, col='green', lty=2)
        
        legend('topright', legend=c('CTA/LPO','TUR/TJE', 'SEA/CA'), col=c('salmon1','slategray3','green'),lty=1)

        abline(v=1.0, col='red', lty=2, lwd=2)

        plot.conf.int(boot.value.dtt, density.dtt, col='slategray3')
        plot.conf.int(boot.value.dcl, density.dcl, col='salmon1')
        plot.conf.int(boot.value.dsc, density.dsc, col='green')

        ## for homozygous derived
        real.value.dcl <- tmp.dcl$hom.der[tmp.dcl$idx == 0]
        real.value.dtt <- tmp.dtt$hom.der[tmp.dtt$idx == 0]
        real.value.dsc <- tmp.dsc$hom.der[tmp.dsc$idx == 0]
        boot.value.dcl <- tmp.dcl$hom.der[tmp.dcl$idx != 0]
        boot.value.dtt <- tmp.dtt$hom.der[tmp.dtt$idx != 0]
        boot.value.dsc <- tmp.dsc$hom.der[tmp.dsc$idx != 0]

        density.dcl <- density(boot.value.dcl)
        density.dtt <- density(boot.value.dtt)
        density.dsc <- density(boot.value.dsc)
        xlims <- range(c(density.dcl$x, density.dtt$x, density.dsc$x))
        ylims <- range(c(density.dcl$y, density.dtt$y, density.dsc$y))

        plot(0,0,type='n', xlim=xlims, ylim=ylims, main=paste(score.class, 'Homozygous derived'), ylab='Density', xlab='Ratio')
        lines(density.dcl, col='salmon1')
        abline(v=real.value.dcl, col='salmon1', lty=2)
        lines(density.dtt, col='slategray3')
        abline(v=real.value.dtt, col='slategray3', lty=2)
        lines(density.dsc, col='green')
        abline(v=real.value.dsc, col='green', lty=2)

        legend('topright', legend=c('CTA/LPO','TUR/TJE','SEA/CA'), col=c('salmon1','slategray3','green'),lty=1)
        abline(v=1.0, col='red', lty=2, lwd=2)
        
        plot.conf.int(boot.value.dtt, density.dtt, col='slategray3')
        plot.conf.int(boot.value.dcl, density.dcl, col='salmon1')
        plot.conf.int(boot.value.dsc, density.dsc, col='green')


    }
}
dev.off()
