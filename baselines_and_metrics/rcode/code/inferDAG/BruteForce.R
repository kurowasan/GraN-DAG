BruteForce <- function(X, scoreName, pars = list(), output = FALSE)
    # Copyright (c) 2011 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{   
    computeStatisticsFromData <- get(paste("computeStatisticsFromData",scoreName,sep=""))
    pars <- computeStatisticsFromData(X,pars)
    p <- dim(X)[2]
    
    if(p > 5)
    {
        stop("This method is intractable for p > 5. p = 6 could probably be done but is not implemented yet.")
    }
    load(paste("../../util_DAGs/allDagsWith",p,"Nodes.RData",sep = ""))
    
    numDags <- dim(allDags)[1]
    if(output)
    {
        cat("In total there are ", numDags, " DAGs\n")        
    }
    allScores <- rep(NA,numDags)
    allScoresStat <- rep(NA,numDags)
    bestScore <- Inf
    bestScoreStat <- Inf
    
    for(i in 1:numDags)
    {
        if(output)
        {
            cat("DAG Number:", i, "\r")
        }
        Adj <- matrix(allDags[i,],p,p)
        initializeState <- get(paste("initializeState",scoreName,sep=""))
        State <- initializeState(X,Adj,pars)
        
        allScores[i] <- State$Score
        allScoresStat[i] <- State$ScoreStat
        if(State$Score < bestScore)
        {
            bestState <- State
            bestScore <- State$Score    
        }
        if(State$ScoreStat < bestScoreStat)
        {
            bestStateStat <- State
            bestScoreStat <- State$ScoreStat    
        }
    }
    
    best <- which.min(allScores)
    
    return(list(bestState = bestState, bestStateStat = bestStateStat, Score = bestScore, allScores = allScores, allScoresStat = allScoresStat, Adj = bestState$Adj))
}

