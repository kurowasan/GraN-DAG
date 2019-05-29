BruteForceFullDags <- function(X, scoreName, pars = list(), output = FALSE)
    # Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{   
    if(scoreName == "SEMGAM")
    {
        aa <- 1
    } else
    {
        computeStatisticsFromData <- get(paste("computeStatisticsFromData",scoreName,sep=""))
        pars <- computeStatisticsFromData(X,pars)
    }
    p <- dim(X)[2]
    
    if(p > 7)
    {
        stop("This method is intractable for p > 7.")
    }
    load(paste("../../util_DAGs/allFullDagsWith",p,"Nodes.RData",sep = ""))
    
    numDags <- dim(allFullDags)[1]
    if(output)
    {
        cat("In total there are ", numDags, " DAGs\n")        
    }
    allScores <- rep(NA,numDags)
    bestScore <- Inf
    
    for(i in 1:numDags)
    {
        if(output)
        {
            cat("DAG Number:", i, "\r")
        }
        Adj <- matrix(allFullDags[i,],p,p)
        initializeState <- get(paste("initializeState",scoreName,sep=""))
        State <- initializeState(X,Adj,pars)
        
        allScores[i] <- State$Score
        if(State$Score < bestScore)
        {
            bestState <- State
            bestScore <- State$Score    
        }
    }
    
    best <- which.min(allScores)
    
    return(list(bestState = bestState, allScores = allScores, Adj = bestState$Adj, score = best))
}
