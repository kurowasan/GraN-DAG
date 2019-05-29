initializeStateSEMGAM <- function(X, initialAdj, pars, penFactor = 1)
    # Copyright (c) 2013  Jonas Peters [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms.
{
    State <- list()
    
    State$n <- dim(X)[1]
    State$p <- dim(X)[2]
    State$Adj <- initialAdj
    State$Res <- matrix(NaN,State$n,State$p)
    State$SingleScores <- rep(0,State$p)
    
    for(node in 1:(State$p))
    {
        parents <- which(State$Adj[,node]==1)
        
        if(length(parents) == 0)
        {
            # use MLE instead of unbiased estimator of the variance
            State$Res[,node] <- X[,node]
        }
        else
        {
            modTmp <- train_gam(X[,parents],X[,node])
            State$Res[,node] <- modTmp$residuals
        }
    }
    
    for(node in 1:(State$p))
    {
        State$SingleScores[node] <- log(var(State$Res[,node]))
    }    
    State$DiffScore <- 0
    State$Score <- sum(State$SingleScores)
    
    return(State)
}
