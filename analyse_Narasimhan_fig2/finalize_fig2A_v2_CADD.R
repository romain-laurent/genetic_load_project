
d <- read.table('counts_fig2A_v2_CADD.txt', header=TRUE)

roh.classes <- c('N','A','B','C')
d$roh.class <- factor(d$roh.class, levels=roh.classes, ordered=TRUE)

pops <- unique(d$pop)
cols <- c('tomato1', 'gold', 'seagreen3', 'lightskyblue1')
max.nb.indiv <- max(table(d$pop)/12)

pdf('fig2A_v2_CADD.pdf', 15, 10, title='Figure 2A v2')

## first plot
neutral <- 'CADD_inf_10'
non.neutral <- 'CADD_sup_28'

res <- matrix(NA, max.nb.indiv, length(pops)*length(roh.classes))
colnames(res) <- c(outer(roh.classes, pops, paste))

for (pop in pops){
    indivs <- unique(d$indiv[d$pop==pop])
    for (idx in 1:length(indivs)){
        tmp.indiv <- indivs[idx]
        for (roh.class in roh.classes){
            tmp.d <- d[d$indiv==tmp.indiv & d$roh.class==roh.class,]
            ratio.neutral <- (2*tmp.d$hom[tmp.d$mut.class==neutral] + tmp.d$het[tmp.d$mut.class==neutral])# / (2*tmp.d$tot[tmp.d$mut.class==neutral])
            ratio.non.neutral <- (2*tmp.d$hom[tmp.d$mut.class==non.neutral] + tmp.d$het[tmp.d$mut.class==non.neutral])# / (2*tmp.d$tot[tmp.d$mut.class==non.neutral])

            tmp <- paste(roh.class, pop)
            res[idx, tmp] <- ratio.non.neutral / ratio.neutral

        }
    }
}

## plot
boxplot.data <- boxplot(res, col=rep(cols, each=length(roh.classes)), ylab=paste(non.neutral,neutral,sep='/'))

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


neutral <- 'CADD_inf_10'
non.neutral <- 'CADD_15_28'

res <- matrix(NA, max.nb.indiv, length(pops)*length(roh.classes))
colnames(res) <- c(outer(roh.classes, pops, paste))

for (pop in pops){
    indivs <- unique(d$indiv[d$pop==pop])
    for (idx in 1:length(indivs)){
        tmp.indiv <- indivs[idx]
        for (roh.class in roh.classes){
            tmp.d <- d[d$indiv==tmp.indiv & d$roh.class==roh.class,]
            ratio.neutral <- (2*tmp.d$hom[tmp.d$mut.class==neutral] + tmp.d$het[tmp.d$mut.class==neutral])# / (2*tmp.d$tot[tmp.d$mut.class==neutral])
            ratio.non.neutral <- (2*tmp.d$hom[tmp.d$mut.class==non.neutral] + tmp.d$het[tmp.d$mut.class==non.neutral])# / (2*tmp.d$tot[tmp.d$mut.class==non.neutral])
            tmp <- paste(roh.class, pop)
            res[idx, tmp] <- ratio.non.neutral / ratio.neutral

        }
    }
}

## plot
boxplot.data <- boxplot(res, col=rep(cols, each=length(roh.classes)), ylab=paste(non.neutral,neutral,sep='/'))

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


neutral <- 'CADD_inf_10'
non.neutral <- 'CADD_10_15'

res <- matrix(NA, max.nb.indiv, length(pops)*length(roh.classes))
colnames(res) <- c(outer(roh.classes, pops, paste))

for (pop in pops){
    indivs <- unique(d$indiv[d$pop==pop])
    for (idx in 1:length(indivs)){
        tmp.indiv <- indivs[idx]
        for (roh.class in roh.classes){
            tmp.d <- d[d$indiv==tmp.indiv & d$roh.class==roh.class,]
            ratio.neutral <- (2*tmp.d$hom[tmp.d$mut.class==neutral] + tmp.d$het[tmp.d$mut.class==neutral])# / (2*tmp.d$tot[tmp.d$mut.class==neutral])
            ratio.non.neutral <- (2*tmp.d$hom[tmp.d$mut.class==non.neutral] + tmp.d$het[tmp.d$mut.class==non.neutral])# / (2*tmp.d$tot[tmp.d$mut.class==non.neutral])
            tmp <- paste(roh.class, pop)
            res[idx, tmp] <- ratio.non.neutral / ratio.neutral

        }
    }
}

## plot
boxplot.data <- boxplot(res, col=rep(cols, each=length(roh.classes)), ylab=paste(non.neutral,neutral,sep='/'))

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
