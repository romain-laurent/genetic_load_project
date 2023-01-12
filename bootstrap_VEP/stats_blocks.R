
if (!('d' %in% ls())){
    d <- read.table('bootstrap_blocks_definition.txt', header=TRUE)
}


par(mfrow=c(1,2))
snps.per.block <- table(d$idx)
plot(density(snps.per.block), main='', xlab='# SNP per block')
#abline(v=quantile(snps.per.block, c(.025,0.975)), col='red')

lengths <- c()
for (i in 1:1000){
    tmp <- d$pos[d$idx==i]
    lengths <- c(lengths, max(tmp)-min(tmp))
}
plot(density(lengths), main='', xlab='Block length (bp)')
