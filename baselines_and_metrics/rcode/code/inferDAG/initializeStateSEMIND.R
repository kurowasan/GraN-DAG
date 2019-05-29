initializeStateSEMIND <- function(X, initialAdj, pars, penFactor = 1)
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
            modTmp <- pars$regr.method(X[,parents], X[,node], pars$regr.pars)
            State$Res[,node] <- modTmp$residuals
        }
    }
    
    for(node in 1:(State$p))
    {
        indtestres <- pars$indtest.method(State$Res[,-node],State$Res[,node])
        State$SingleScoresStat[node] <- indtestres$statistic        
        State$SingleScores[node] <- - max(log(indtestres$p.value)/log(10),-350)
        #State$SingleScores[node] <- - log(indtestres$p.value)/log(10)
    }    
    State$DiffScore <- 0
    cS <- computeScoreSEMIND(State,pars,penFactor)
    State$Score <- cS$score
    State$ScoreStat <- cS$scoreStat
    
    return(State)
}
