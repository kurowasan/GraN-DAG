computeScoreCAMGauss <- function(X, Adj, fullDAG = TRUE)
    # Copyright (c) 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    if(fullDAG)
    {
        a <- computeCausOrder(G)
        for(i in 1:(p-1))
        {
            G[a[i],a[(i+1):p]] <- rep(1,p-i)
        }
    }
    
    
    eachResVar <- rep(0,p)
    
    
    for(node in 1:(p))
    {
        parents <- which(Adj[,node]==1)
        
        if(length(parents) == 0)
        {
            # use MLE instead of unbiased estimator of the variance
            eachResVar[node] <- var(X[,node])
        }
        else
        {
            mod <- train_gam(X[,parents], X[,node])
            # use MLE instead of unbiased estimator of the variance
            eachResVar[node] <- var(mod$residuals)
        }
    }
    
    # the following is the negloglikelihood
    # negloglik <- n/2 * log((2*pi)^p*prod(eachResVar)) + n/2 * p
    
    # for comparison we compute sth else
    loglik <- -sum(log(eachResVar))
    
    # minimize this BIC score!
    #score <- negloglik + penalization
    return(loglik)
}
