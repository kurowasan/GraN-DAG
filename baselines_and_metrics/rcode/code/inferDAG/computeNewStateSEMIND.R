computeNewStateSEMIND <- function(oldState,index,X,pars,penFactor,perm = 0, output = FALSE)
    # Copyright (c) 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    n <- dim(X)[1]
    p <- dim(X)[2]
    i <- index[1]
    j <- index[2]
    
    # take over oldState and change Adj matrix
    newState <- oldState
    newState$Adj <- oldState$Adj
    
    if(perm==0)
    {
        newState$Adj[i,j] <- (newState$Adj[i,j] + 1) %% 2
        
        if(newState$Adj[j,i] == 1)
        {
            newState$Adj[j,i] <- 0
        #    recomputeNodes <- c(i,j)
            recomputeNodes <- c(i,j)
        }
        else
        {
        #    recomputeNodes <- c(j)
            recomputeNodes <- 1:p
        }
    }
    else
    {
        newState$CausOrder <- oldState$CausOrder
    }
    if(perm==1)
    {
        # exchange i and j
        tmpp <- 1:p
        tmpp[c(newState$CausOrder[i],newState$CausOrder[j])] <- tmpp[c(newState$CausOrder[j],newState$CausOrder[i])]
        newState$Adj[tmpp,tmpp] <- newState$Adj
        newState$CausOrder[c(i,j)] <- newState$CausOrder[c(j,i)]
        recomputeNodes <- newState$CausOrder[i:j]
        
    }
    if(perm ==2)
    {
        # put i to position j 
        tmpi <- newState$CausOrder[i]
        if(i>j)
        {
            newState$CausOrder[(j+1):i] <- newState$CausOrder[j:(i-1)]
        }
        else
        {
            newState$CausOrder[i:(j-1)] <- newState$CausOrder[(i+1):j]
        }
        newState$CausOrder[j] <- tmpi
        newState$Adj <- upper.tri(matrix(-1, p, p), diag = 0)
        newState$Adj <- newState$Adj[newState$CausOrder,newState$CausOrder]
        
        recomputeNodes <- newState$CausOrder[i:j]
    }
    
    
    for(node in recomputeNodes)
    {
        parents <- which(newState$Adj[,node]==1)
        
        if(length(parents) == 0)
        {
            newState$Res[,node] <- X[,node]
        }
        else
        {
            modTmp <- pars$regr.method(X[,parents], X[,node], pars$regr.pars)
            newState$Res[,node] <- modTmp$residuals
        }
    }
    for(node in recomputeNodes)
    {
        indtestres <- pars$indtest.method(newState$Res[,-node],newState$Res[,node])
        newState$SingleScoresStat[node] <- indtestres$statistic
        newState$SingleScores[node] <- - max(log(indtestres$p.value)/log(10),-350)
        #newState$SingleScores[node] <- - log(indtestres$p.value)/log(10)
    }    
    
    diffScore <- 0
    #    for(node in c(i,j))
    #    {
    #        indtestres <- pars$indtest.method(newState$Res[,-node],newState$Res[,node])
    #        newState$SingleScores[node] <- - indtestres$p.value
    #        diffScore <- diffScore + newState$SingleScores[node] - oldState$SingleScores[node]
    #    }
    newState$DiffScore <- diffScore
    
    cS <- computeScoreSEMIND(newState, pars, penFactor)
    newState$Score <- cS$score
    newState$ScoreStat <- cS$scoreStat 
    
    #
    # problem with maximum score: you do not realize if you improve on one node that is not
    # dealing with the most dependent residuals.
    
    return(newState)
}
