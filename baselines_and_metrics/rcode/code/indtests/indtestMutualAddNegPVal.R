indtestMutualAddNegPVal <- function(X,indtest,parss = list())
{
    dimension <- dim(X)[2]
    pValVec <- rep(1,dimension-1)
    statVec <- rep(-3,dimension-1)
    
    j <- 1
    for(i in (dimension:2))
    {
        resTmp <- indtest(X[,1:(i-1)],X[,i],pars = parss)
        pValVec[j] <- resTmp$p.value
        statVec[j] <- resTmp$statistic        
        j <- j + 1
    }       
    ind <- which.min(pValVec)
    stat <- statVec[ind]
    # bonferroni
    pval = pValVec[ind] / (j-1)
    res <- list(statistic = stat, p.value = pval)
    return(res)
}
