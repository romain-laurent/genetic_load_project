

d <- read.table('counts_fig2A.txt', header=TRUE)
roh.len <- read.table('ROH_lengths.txt', header=TRUE)
annotations <- read.table('../VEP/VEP_annotations.txt', header=TRUE)

syn <- 'VEP_SS'
lof <- 'VEP_LOF'

syn.A <- sprintf('%s_C', syn)
syn.NA <- sprintf('%s_N', syn)
lof.A <- sprintf('%s_C', lof)
lof.NA <- sprintf('%s_N', lof)
    
global.counts <- table(annotations$class)
global.ratio <- global.counts[syn] / global.counts[lof]


pops <- unique(d$pop)
cols <- c('tomato1', 'gold', 'seagreen3', 'lightskyblue1')

res <- NULL
pop.max <- max(table(d$pop))
pvals <- c()
for (pop in pops){
    tmp.d <- d[d$pop==pop,]
    NA.ratio <- tmp.d[,syn.NA] / tmp.d[,lof.NA]
    A.ratio <- tmp.d[,syn.A] / tmp.d[,lof.A]
    NA.ratio <- c(NA.ratio, rep(NA, pop.max - length(NA.ratio)))
    A.ratio <- c(A.ratio, rep(NA, pop.max - length(A.ratio)))
    res[[pop]] <- data.frame(non.autozygous=NA.ratio, autozygous=A.ratio)
    cat(pop,'\n')
#    print(wilcox.test(NA.ratio, A.ratio, paired=TRUE))
    pvals <- c(pvals, wilcox.test(NA.ratio, A.ratio, paired=TRUE)$p.value)
}

pdf('fig2A_v1.pdf', 18, 9, title='Ratio Syn/LOF')
boxplot(unlist(res, recursive=FALSE), col=rep(cols, each=2), ylab='Syn/LOF ratio' )
abline(h=global.ratio, col='red', lty=2)

y.max <- 80
y.min <- y.max - 2
for (i in 1:4){
    x.max <- 2*i
    x.min <- x.max - 1
    lines(c(x.min, x.min, x.max, x.max), c(y.min, y.max, y.max, y.min))
    text(x.min+.5, y.max, sprintf('p = %.2e', pvals[i]), pos=3)
}

## NA.ratio <- d[,syn.NA] / d[,lof.NA]
## A.ratio <- d[,syn.A] / d[,lof.A]
## tmp <- data.frame(non.autozygous=NA.ratio, autozygous=A.ratio)
## boxplot(tmp, ylab='Syn/LOF ratio')
## abline(h=global.ratio, col='red', lty=2)
