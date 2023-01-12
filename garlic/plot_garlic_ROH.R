
pops <- c('CTA','LPO','TUR','TJE')
roh.classes <- c('A','B','C')
win.size <- 50


## read data once
if (!('final.d' %in% ls())){
    final.d <- NULL
    for (roh.class in roh.classes){
        
        d <- read.table(sprintf('ROH_class_%s_winsize_%d.txt', roh.class, win.size), header=TRUE)
        for (pop in pops){
            indivs <- unique(d$indiv[d$pop==pop])
            for (indiv in indivs){
                tmp.d <- d[d$indiv == indiv,]
                size <- sum(tmp.d$stop - tmp.d$start)
                line <- data.frame(pop=pop, indiv=indiv, roh.class=roh.class, size=size)
                final.d <- rbind(final.d, line)
            }
        }
    }
    final.d$size <- final.d$size / 1e6
}

cols <- c('tomato1', 'gold', 'seagreen3', 'lightskyblue1')
names(cols) <- pops

make.transparent <- function(colname, alphas){
    alphas <- alphas * 255 / 100
    tmp <- col2rgb(colname)
    new.cols <- c()
    for (alpha in alphas){
        new.cols <- c(new.cols, rgb(red=tmp[1], green=tmp[2], blue=tmp[3],alpha=alpha, maxColorValue=255))
    }
    new.cols
}

pdf('ROH_per_indiv.pdf', 15, 9, title='ROH GARLIC per indiv')
par(mfrow=c(2,2))
ymax <- max(final.d$size)
for (pop in pops){
    tmp.d <- final.d[final.d$pop==pop,]
    tmp.col <- cols[pop]
    my.cols <- make.transparent(tmp.col, alphas=c(100,50,20))
    barplot(size~roh.class+indiv, data=tmp.d, beside=TRUE, ylim=c(0,ymax), main=pop, xlab='', ylab='Total ROH length (Mb)', col=my.cols, las=2)
    legend('topleft', legend=c('A','B','C'), fill=my.cols, horiz=TRUE)
}
dev.off()


pdf('correlation_HBD_ROH.pdf', 15, 10, title='Correlation GARLIC/FEstim')
par(mfrow=c(2,2))
for (pop in pops){
    d.hbd <- read.table(sprintf('../Festim/results/flod-%s/HBD.segments.txt', pop), header=TRUE)
    d.hbd$len.bp <- (d.hbd$FIN_bp - d.hbd$DEB_bp) / 1e6
    tmp.d <- NULL
    for (indiv in unique(final.d$indiv[final.d$pop==pop])){
        total.hbd <- sum(d.hbd$len.bp[d.hbd$IID==indiv])
        total.roh <- final.d$size[final.d$indiv==indiv & final.d$roh.class == 'C']
        tmp.d <- rbind(tmp.d, data.frame(indiv=indiv, hbd=total.hbd, roh=total.roh))
    }


    model <- lm(roh~hbd, data=tmp.d)
    intercept <- model$coefficients[1]
    slope <- model$coefficients[2]

    main <- sprintf('%s -- ROH=%.2f HBD + %.2f', pop, slope, intercept)
    plot(tmp.d$hbd, tmp.d$roh, main=main, xlab='Total length of HBD segments (Mb)', ylab='Total length of ROH class C segments (Mb)')
    abline(model, col='red', lty=2)


}
dev.off()
