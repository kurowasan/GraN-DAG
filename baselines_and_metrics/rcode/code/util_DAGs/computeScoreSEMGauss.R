computeScoreSEMGauss <- function(X, Adj)
    # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    SigmaHat <- cov(X)
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    I <- diag(rep(1,p))
    B <- diag(rep(0,p))
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
            mod <- lm(X[,node] ~ X[,parents])  
            # use MLE instead of unbiased estimator of the variance
            eachResVar[node] <- var(mod$residuals)
            B[node,parents] <- mod$coef[2:(length(parents) + 1)]
        }
    }
    
    numPars <- (sum(B!=0) + p)
    
    # sigmaHatSq <- 1/(pn) sum_k sum(State$eachResVar[k] - overallMean)^2
    # but all the means are by cosntruction of the regression zero.
    SigmaNInvHat <- diag(1/eachResVar)
     
    # the following is the negloglikelihood
    #negloglik <- n/2 * log((2*pi)^p*prod(eachResVar)) + n/2 * sum( diag( t(I-B) %*% (I-B) %*% SigmaNInvHat %*% SigmaHat) )
    negloglik <- n/2 * log((2*pi)^p*prod(eachResVar)) + n/2 * p
    penalization <- log(n)/2 * numPars
   
    if(1==0)
    {
        show("eachResVar")
        show(eachResVar)
        show("negloglik SEM")
        show(negloglik)
        show("pen SEM")
        show(penalization)
    }    
    # minimize this BIC score!
    score <- negloglik + penalization
    return(score)
}
