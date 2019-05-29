initializeStateSEMSEV <- function(X, initialAdj, pars, penFactor = 1)
# Copyright (c) 2012-2012  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.
{
    State <- list()

    State$n <- dim(X)[1]
    State$p <- dim(X)[2]
    State$Adj <- initialAdj
    State$B <- diag(rep(0,State$p))

    for(node in 1:(State$p))
    {
        parents <- which(State$Adj[,node]==1)
        
        if(length(parents) == 0)
        {
            # use MLE instead of unbiased estimator of the variance
            State$eachResVar[node] <- (State$n-1)/State$n*var(X[,node])
        }
        else
        {
            mod <- lm(X[,node] ~ X[,parents])  
            # use MLE instead of unbiased estimator of the variance
            State$eachResVar[node] <- (State$n-1)/State$n*var(mod$residuals)
            State$B[node,parents] <- mod$coef[2:(length(parents) + 1)]
        }
    }
    
    State$sumResVar <- sum(State$eachResVar)
    # number of pars = number of non-zero coefs + 1 (noise variance)
    State$numPars <- (sum(State$B!=0) + 1)
    State$Score <- computeScoreSEMSEV(State,pars,penFactor)
    
    return(State)
}
