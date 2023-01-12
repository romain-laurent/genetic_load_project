
pops <- c('CTA','LPO','TUR','TJE')
win.size <- 1000
nb.line <- 0
for (pop in pops){
    filename <- sprintf('diversity_%s_win%dk.windowed.pi', pop, win.size)
    command <- sprintf('wc -l %s', filename)
    res <- system(command, intern=TRUE)
    res <- as.integer(strsplit(res, ' ')[[1]][1])
    nb.line <- max(nb.line, res)
}

divs <- matrix(NA, nrow=length(pops), ncol=nb.line)
rownames(divs) <- pops

for (pop in pops){
    d <- read.table(sprintf('diversity_%s_win%dk.windowed.pi', pop, win.size), header=TRUE)
    divs[pop,1:nrow(d)] <- d$PI
}

cols <- c('tomato1', 'gold', 'seagreen3', 'lightskyblue1')
names(cols) <- pops
pdf('nucleotide_diversity.pdf', 10, 7, title='Nucleotide diversity')
boxplot(t(divs), col=cols, ylab=sprintf('Nucleotide diversity computed in %d kbp regions',win.size))

