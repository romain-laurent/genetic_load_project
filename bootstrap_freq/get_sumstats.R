
pops <- c('CTA','LPO','TUR','TJE')

if (!('d' %in% ls())){
    d <- read.table('precomputed_freq_stats.txt', header=TRUE)
    pop.table <- NULL
    for (pop in pops){
        tmp <- read.table(sprintf('../samples_%s.txt',pop))
        colnames(tmp) <- c('indiv')
        tmp$pop <- pop
        pop.table <- rbind(pop.table, tmp)
    }
}

d$pop <- pop.table$pop[match(d$indiv, pop.table$indiv)]
scores <- c("CADD_inf_10", "CADD_10_15", "CADD_15_28", "CADD_sup_28", "total")

res <- data.frame(pop=rep(pops, each=length(scores)), score=rep(scores, length(pops)), mean.der=NA, sd.der=NA, mean.hom.der=NA, sd.hom.der=NA)
scores <- scores[1:(length(scores)-1)]

for (pop in pops){
    indivs <- pop.table$indiv[pop.table$pop==pop]
    total.derived <- rep(0, length(indivs))
    total.hom.derived <- rep(0, length(indivs))
    
    for (score in scores){
        wanted.derived <- c()
        wanted.hom.derived <- c()
        for (indiv in indivs){
            tmp <- d[d$indiv==indiv & d$score.class==score,]
            wanted.derived <- c(wanted.derived, sum(tmp$het + 2*tmp$hom.der))
            wanted.hom.derived <- c(wanted.hom.derived, sum(tmp$hom.der))
        }
        total.derived <-  total.derived + wanted.derived
        total.hom.derived <- total.hom.derived + wanted.hom.derived
        
        row <- which(res$pop==pop & res$score==score)
        res$mean.der[row] <- mean(wanted.derived)
        res$sd.der[row] <- sd(wanted.derived)
        res$mean.hom.der[row] <- mean(wanted.hom.derived)
        res$sd.hom.der[row] <- sd(wanted.hom.derived)
    }
    score <- 'total'
    row <- which(res$pop==pop & res$score==score)
    res$mean.der[row] <- mean(total.derived)
    res$sd.der[row] <- sd(total.derived)
    res$mean.hom.der[row] <- mean(total.hom.derived)
    res$sd.hom.der[row] <- sd(total.hom.derived)
    
}


save(res, file='sumstats.Rdata')
