
## read data
rab <- read.table('Narasimhan_table_S6.csv', header=TRUE, sep='\t', row.names=1)
inbreeding <- read.table('Gazal_supp.csv', header=TRUE, sep='\t')



## modify rab table
## remove BB row
rab <- rab[rownames(rab)!='BB',]
## add empty YRI row
rab['YRI',] <- NA
## compute lower triangular values
lower <- 1/t(rab)
## put them in place
rab[lower.tri(rab)] <- lower[lower.tri(lower)]


## compute sumstat for inbreeding in each pop
pops <- colnames(rab)
final.inbreeding <- rep(-1, length(pops))
names(final.inbreeding) <- pops
for (pop in pops){
    tmp.d <- inbreeding$f[inbreeding$POP == pop]
    final.inbreeding[pop] <- length(which(tmp.d > 0)) / length(tmp.d)
}

cols <- rep('orange', length(pops))
tmp.cols <- c('red','green','blue','purple')
names(tmp.cols) <- c('AFR','EUR','EAS','SAS')
names(cols) <- pops
for (pop in pops){
    tmp <- unique(inbreeding$SUPER.POP[inbreeding$POP==pop])[1]
    cols[pop] <- tmp.cols[tmp]
}

## plot
pdf('correlation.pdf', 15, 15, title='Correlation')
par(mfrow=c(4,5))
for (ref.pop in pops){
    print(ref.pop)
    xlims <- c(min(final.inbreeding), max(final.inbreeding)*1.2)
    plot(final.inbreeding, rab[ref.pop,], main=ref.pop, xlab='Mean f', ylab='Rab',xlim=xlims )
    print(cor.test(c(final.inbreeding), as.numeric(rab[ref.pop,]), method='pearson', na.rm=T))
   text(final.inbreeding, rab[ref.pop,], names(final.inbreeding),pos=4)
}
dev.off()



pdf('correlation_ref_CEU.pdf', 7, 7, title='Correlation Rab / inbreeding')
ref.pop <- 'CEU'
ylims <- c(min(rab[ref.pop,], na.rm=TRUE), max(rab[ref.pop,], na.rm=TRUE)*1.01)
plot(final.inbreeding, rab[ref.pop,], xlab='Mean f', ylab='Rab',ylim=ylims, col=cols[names(final.inbreeding)])
text(final.inbreeding, rab[ref.pop,], names(final.inbreeding), pos=3, col=cols[names(final.inbreeding)])

dev.off()
