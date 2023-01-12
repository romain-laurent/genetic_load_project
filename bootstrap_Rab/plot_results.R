
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


dcl <- read.table('bootstrap_Rab_CTA_LPO.txt', header=TRUE)
dtt <- read.table('bootstrap_Rab_TUR_TJE.txt', header=TRUE)
dsc <- read.table('bootstrap_Rab_SEA_CA.txt', header=TRUE)

pdf('bootstrap_Rab.pdf', 15, 10, title='Bootstrap Rab')
for (i in 1:length(all.scores)){
    scores <- all.scores[[i]]
    par(mfcol=c(length(scores), 2))
    ## plot raw stat (original Rxy from Do)
    for (score.class in scores){
        
        real.value.dcl <- dcl[dcl$idx == 0, score.class]
        real.value.dtt <- dtt[dtt$idx == 0, score.class]
        real.value.dsc <- dsc[dsc$idx == 0, score.class]
        boot.value.dcl <- dcl[dcl$idx != 0, score.class]
        boot.value.dtt <- dtt[dtt$idx != 0, score.class]
        boot.value.dsc <- dsc[dsc$idx != 0, score.class]

        density.dcl <- density(boot.value.dcl)
        density.dtt <- density(boot.value.dtt)
        density.dsc <- density(boot.value.dsc)
        xlims <- range(c(density.dcl$x, density.dtt$x, density.dsc$x))
        ylims <- range(c(density.dcl$y, density.dtt$y, density.dsc$y))

        plot(0,0,type='n', xlim=xlims, ylim=ylims, main=score.class, ylab='Density', xlab='Ratio')
        lines(density.dcl, col='salmon1')
        abline(v=real.value.dcl, col='salmon1', lty=2)
        lines(density.dtt, col='slategray3')
        abline(v=real.value.dtt, col='slategray3', lty=2)
    #    lines(density.dsc, col='green')
   #     abline(v=real.value.dsc, col='green', lty=2)
        
        legend('topright', legend=c('CTA/LPO','TUR/TJE','SEA/CA'), col=c('salmon1','slategray3','green'),lty=1)

        abline(v=1.0, col='red',lty=2, lwd=2)
        
        plot.conf.int(boot.value.dtt, density.dtt, col='slategray3')
        plot.conf.int(boot.value.dcl, density.dcl, col='salmon1')
  #      plot.conf.int(boot.value.dsc, density.dsc, col='green')

    }

    ## plot Rab stat from Narasimhan, i.e. R'xy from Do
    for (score.class in scores){
        neutral.class <- levels(scores)[1]
        real.value.dcl <- dcl[dcl$idx == 0, score.class] / dcl[dcl$idx == 0, neutral.class]
        real.value.dtt <- dtt[dtt$idx == 0, score.class] / dtt[dtt$idx == 0, neutral.class]
        real.value.dsc <- dsc[dsc$idx == 0, score.class] / dsc[dsc$idx == 0, neutral.class]
        boot.value.dcl <- dcl[dcl$idx != 0, score.class] / dcl[dcl$idx != 0, neutral.class]
        boot.value.dtt <- dtt[dtt$idx != 0, score.class] / dtt[dtt$idx != 0, neutral.class]
        boot.value.dsc <- dsc[dsc$idx != 0, score.class] / dsc[dsc$idx != 0, neutral.class]

        density.dcl <- density(boot.value.dcl)
        density.dtt <- density(boot.value.dtt)
        density.dsc <- density(boot.value.dsc)
        xlims <- range(c(density.dcl$x, density.dtt$x, density.dsc$x))
        ylims <- range(c(density.dcl$y, density.dtt$y, density.dsc$y))

        plot(0,0,type='n', xlim=xlims, ylim=ylims, main=score.class, ylab='Density', xlab='Ratio')
        lines(density.dcl, col='salmon1')
        abline(v=real.value.dcl, col='salmon1', lty=2)
        lines(density.dtt, col='slategray3')
        abline(v=real.value.dtt, col='slategray3', lty=2)
#        lines(density.dsc, col='green')
#        abline(v=real.value.dsc, col='green', lty=2)

        legend('topright', legend=c('CTA/LPO','TUR/TJE','SEA/CA'), col=c('salmon1','slategray3','green'),lty=1)
        abline(v=1.0, col='red',lty=2, lwd=2)
        
        plot.conf.int(boot.value.dtt, density.dtt, col='slategray3')
        plot.conf.int(boot.value.dcl, density.dcl, col='salmon1')
 #       plot.conf.int(boot.value.dsc, density.dsc, col='green')
    }

}
dev.off()
