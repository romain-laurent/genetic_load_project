
pops <- c('CTA','LPO','TUR','TJE')
roh.classes <- c('A','B','C')

nb.group <- 2


## read data once
if (!('final.d' %in% ls())){
    final.d <- NULL
    d <- read.table(sprintf('ROH_classes_%d_groups.txt',nb.group),header=TRUE)
    for (pop in pops){
        indivs <- unique(d$indiv[d$pop==pop])
        for (roh.class in roh.classes){
            for (indiv in indivs){
                tmp.d <- d[d$indiv == indiv & d$class==roh.class,]
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

pdf(sprintf('ROH_per_indiv_reclassed_%d_groups.pdf', nb.group), 15, 9, title=sprintf('ROH GARLIC per indiv reclassed %d groups', nb.group))
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


pdf(sprintf('correlation_HBD_ROH_reclassed_%d_groups.pdf', nb.group), 15, 10, title=sprintf('Correlation GARLIC/FEstim reclassed %d groups', nb.group))
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
