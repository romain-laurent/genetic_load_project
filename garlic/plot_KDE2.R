


pops <- c('CTA','LPO','TUR','TJE')
thresholds <- NULL
winsizes <- seq(30,100,by=10)
values <- c()

pdf('KDEs_force_winsize.pdf', 20, 80)
par(mfcol=c(length(winsizes),4))
for (pop in pops){
    for (winsize in winsizes){
        filename <- sprintf('results_force_winsize/ROH_%s_winsize_%d.%dSNPs.kde', pop, winsize, winsize)
        d <- read.table(filename)
        plot(d, type='l', main=sprintf('%s %d', pop, winsize))
        filename <- sprintf('results_force_winsize/ROH_%s_winsize_%d.log', pop, winsize)
        command <- sprintf('grep LOD %s | grep Selected | cut -d" " -f5', filename)
        value <- as.numeric(system(command, intern=T))
        values <- c(values, value)
    }
}
dev.off()


values <- data.frame(pop=rep(pops, each=length(winsizes)),
                     winsize=rep(winsizes, length(pops)),
                     value=values)

pdf('LOD_cutoff_force_winsize.pdf', 10, 10)
par(mfrow=c(1,1))
col <- 1
plot(values$winsize, values$value, type='n', xlab='Window size', ylab='LOD score cutoff')
for (pop in pops){
    tmp <- values[values$pop==pop,]
    lines(tmp$winsize, tmp$value, col=col, lwd=2)
    col <- col+1
}

legend('bottomleft', pops, col=1:length(pops), lwd=2)
dev.off()
