prefixes <- c('CTA','LPO','TUR','TJE')

pdf('SFS.pdf', 15, 10)
for (prefix in prefixes){
    cat(prefix,'\n')
    d <- read.table(sprintf('SFS_%s.frq.count', prefix), skip=1)
    to.keep <- d[,6] != 0 
    d <- d[to.keep,]
    to.plot <- table(d[,6])
    plot(to.plot, xlab='# ALT allele', ylab='Count', main=paste(prefix, nrow(d)))
    lines(1:length(to.plot), to.plot[1] / 1:length(to.plot), col='red', type='b')
}
dev.off()
