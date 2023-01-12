
## read data
d <- read.table('../analyse_Narasimhan_fig2/counts_fig2A.txt', header=TRUE)

## define constants
pops <- unique(d$pop)
cols <- c('tomato1', 'gold', 'seagreen3', 'lightskyblue1')

max.nb.indiv <- max(table(d$pop))
roh.classes <- c('N','A','B','C')

## compute expected ratio
all.annotations <- read.table('../VEP/VEP_annotations.txt', header=TRUE)
all.counts <- table(all.annotations$class)


pdf('boxplot_ratios.pdf', title='Boxplots ratios', 15, 15)

par(mfrow=c(2,1))

## first plot -> non neutral = LoF
neutral <- 'VEP_SS'
not.neutral <- 'VEP_LOF'

global.ratio <- all.counts[not.neutral] / all.counts[neutral]

res <- matrix(NA, max.nb.indiv, length(pops)*length(roh.classes))
colnames(res) <- c(outer(roh.classes, pops, paste))
## compute ratios for each indiv
for (pop in pops){
    tmp.d <- d[d$pop==pop,]
    for (roh.class in roh.classes){
        colname <- paste(roh.class, pop)
        neutral.str <- sprintf('%s_%s', neutral, roh.class)
        not.neutral.str <- sprintf('%s_%s', not.neutral, roh.class)

        ratio <- tmp.d[,not.neutral.str] / tmp.d[,neutral.str]
       #ratio <- tmp.d[,neutral.str] / tmp.d[,not.neutral.str]
        res[1:length(ratio),colname] <- ratio
    }
    
}

## plot
boxplot.data <- boxplot(res, col=rep(cols, each=length(roh.classes)), ylab='LoF/Syn ratio')
abline(h=global.ratio, lty=2, col='red')

## add pvalues
tmp.x <- 0
for (pop in pops){
    non.auto.str <- paste('N', pop, sep=' ')
    non.auto <- res[,non.auto.str]
    tmp.x <- tmp.x + 1
    for (roh.class in roh.classes[2:length(roh.classes)]){
        tmp.x <- tmp.x + 1
        auto.str <- paste(roh.class, pop, sep=' ')
        auto <- res[,auto.str]
        pval <- wilcox.test(non.auto, auto, paired=TRUE)$p.value
        text(tmp.x, boxplot.data$stats[5,tmp.x], sprintf('%.2e', pval), pos=3)
    }
}

######################################################################
## second plot -> non neutral = non syn
neutral <- 'VEP_SS'
not.neutral <- 'VEP_NS'

global.ratio <- all.counts[not.neutral] / all.counts[neutral]

res <- matrix(NA, max.nb.indiv, length(pops)*length(roh.classes))
colnames(res) <- c(outer(roh.classes, pops, paste))
## compute ratios for each indiv
for (pop in pops){
    tmp.d <- d[d$pop==pop,]
    for (roh.class in roh.classes){
        colname <- paste(roh.class, pop)
        neutral.str <- sprintf('%s_%s', neutral, roh.class)
        not.neutral.str <- sprintf('%s_%s', not.neutral, roh.class)

        ratio <- tmp.d[,not.neutral.str] / tmp.d[,neutral.str]
#       ratio <- tmp.d[,neutral.str] / tmp.d[,not.neutral.str]
        res[1:length(ratio),colname] <- ratio
    }
    
}

## plot
boxplot.data <- boxplot(res, col=rep(cols, each=length(roh.classes)), ylab='NonSyn/Syn ratio')
abline(h=global.ratio, lty=2, col='red')

## add pvalues
tmp.x <- 0
for (pop in pops){
    non.auto.str <- paste('N', pop, sep=' ')
    non.auto <- res[,non.auto.str]
    tmp.x <- tmp.x + 1
    for (roh.class in roh.classes[2:length(roh.classes)]){
        tmp.x <- tmp.x + 1
        auto.str <- paste(roh.class, pop, sep=' ')
        auto <- res[,auto.str]
        pval <- wilcox.test(non.auto, auto, paired=TRUE)$p.value
        text(tmp.x, boxplot.data$stats[5,tmp.x], sprintf('%.2e', pval), pos=3)
    }
}

dev.off()
