
## computes GMM and boundaries between ROH classes the same way GARLIC does
## difference in results with GARLIC are caused by machine precision and differences in random initializations of algorithms

library(mclust)
groups <- list(c('CTA','LPO'), c('TUR','TJE'))
#groups <- list(c('CTA','LPO','TUR','TJE'))

## function we will need to solve to find classes boundaries
to.solve <- function(x, pars){
    pars$prop1 * dnorm(x, pars$mean1, pars$std1) - pars$prop2 * dnorm(x, pars$mean2, pars$std2)
}


## read all data once
if (!('d' %in% ls())){
    d <- read.table('all_ROH_all_pops.txt', header=TRUE)
}


final.d <- NULL
## now we work per group
for (idx.group in 1:length(groups)){
    tmp.group <- groups[[idx.group]]
    cat(tmp.group, '\n')

    ## we select the data that correspond to the group we are working on
    group.d <- NULL
    for (pop in tmp.group){
        tmp.d <- d[d$pop == pop,]
        group.d <- rbind(group.d, tmp.d)
    }

    ## we compute the GMM
    fit <- Mclust(group.d$cM, G=3, modelNames='V')
    ## we get the results
    props <- fit$parameters$pro
    means <- fit$parameters$mean
    vars <- fit$parameters$variance$sigmasq
    
    tmp <- data.frame(prop=props, mean=means,var=vars, row.names=c('A','B','C'))
    print(tmp)

#    plot(density(group.d$cM, from=0), xlim=c(0,3))
    
    
    ## we find the boundaries
    bounds <- c()
    for (idx.bound in 2:3){
        params <- list(prop1=props[idx.bound-1],
                       mean1=means[idx.bound-1],
                       std1=sqrt(vars[idx.bound-1]),
                       prop2=props[idx.bound],
                       mean2=means[idx.bound],
                       std2=sqrt(vars[idx.bound]))
        res <- uniroot(to.solve, pars=params, interval=c(params$mean1,params$mean2),tol=1e-8)
        ## just a check that we converged
        if (res$iter > 900){
            cat('uniroot probably did not converge... bounds may be fucked...\n')
        }
        bounds <- c(bounds, res$root)
    }
    cat('boundaries:', bounds,'\n')

    ## we assign ROHs to determined classes
    group.d$class <- rep('0', nrow(group.d))
    group.d$class[group.d$cM < bounds[1]] <- 'A'
    group.d$class[group.d$cM > bounds[2]] <- 'C'
    group.d$class[group.d$class=='0'] <- 'B'
    final.d <- rbind(final.d, group.d)
}

write.table(final.d, file=sprintf('ROH_classes_%d_groups.txt', length(groups)), quote=FALSE, row.names=FALSE)



