oneStepCheckingUntilFirstMinK <- function(StateOld,scoreName,pars,X,k,penFactor,checkDAGs)
    # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    p <- StateOld$p
    if(scoreName == "SEMSEV")
    {
        # start with largest residual variance.
        # sortNodes <- sort(StateOld$eachResVar, index.return=TRUE, decreasing=TRUE)$ix
        varTmp <- StateOld$eachResVar - min(StateOld$eachResVar)
        #varTmp <- abs(StateOld$eachResVar - mean(StateOld$eachResVar))
        varTmp <- varTmp + min(varTmp[varTmp>0])
        varTmp <- varTmp / sum(varTmp)
        sortNodes <- sample(x = 1:p, size = p, p = varTmp, replace = FALSE)
    } 
    if(scoreName == "SEMIND")
    {
        sortNodes <- sample(x = 1:p, size = p, p = pmax(StateOld$SingleScores,0.000001), replace = FALSE)   
    } else
    {
        sortNodes <- 1:p
    }
    
    #cat("\n")
    #show(sortNodes)
    
    # collect all different neighbors
    index <- array(dim = c(p*(p-1),2))
    indexCount <- 0
    for(i in 1:p)
    {
        for(j in 1:(p-1))
        {
            # do not add sth on the diagonal
            j <- j + (j>(i-1))
            indexCount <- indexCount + 1
            index[indexCount,] <- c(sortNodes[j],sortNodes[i])
        }
    }
    
    # permute this list randomly
    #    index <- index[sample(p*(p-1), replace = FALSE),]
    
    TempState <- StateOld
    TempState$Score <- Inf
    TempState$ScoreStat <- Inf
    
    State1 <- StateOld
    whichStep1 <- c(-1, -1)
    madeStep1 <- FALSE
    indexCount <- 0
    tried <- 0
    while( ((madeStep1 == FALSE)|(tried < k)) & (indexCount < (p*(p-1))) )
    {    
        indexCount <- indexCount + 1
        i <- index[indexCount,1]
        j <- index[indexCount,2]
        
        candidateAdj <- StateOld$Adj
        #        # no reversals
        candidateAdj[i,j] <- (candidateAdj[i,j] + 1) %% 2
        candidateAdj[j,i] <- 0
        
        # no cycles
        if( (!containsCycle(candidateAdj)) & ((i != StateOld$tabuInd[1]) | (j != StateOld$tabuInd[2])) ) 
        {
            tried <- tried +1
            checkDAGs <- checkDAGs + 1
            computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
            newState <- computeNewState(StateOld,c(i,j),X,pars,penFactor)
            
            if( (scoreName == "SEMIND") && (max(StateOld$SingleScores) == 350) )
            {
                if(newState$ScoreStat < State1$ScoreStat)
                {
                    State1 <- newState
                    madeStep1 <- TRUE
                    State1$tabuInd <- c(i,j)
                    whichStep1 <- c(i,j)
                }
                if(newState$ScoreStat < TempState$ScoreStat)
                {
                    TempState <- newState
                    TempState$tabuInd <- c(i,j)
                }
            } else
            {
                if(newState$Score < (1.0 * State1$Score))
                {
                    State1 <- newState
                    madeStep1 <- TRUE
                    State1$tabuInd <- c(i,j)
                    whichStep1 <- c(i,j)
                }
                if(newState$Score < TempState$Score)
                {
                    TempState <- newState
                    TempState$tabuInd <- c(i,j)
                }
            }
        }
    }
    if((scoreName == "SEMIND") && (0 == 1))
    {
        cat("\n- log(p-values) of the different nodes:\n")
        show(State1$SingleScores)
    }
    
    return(list(State = State1, madeStep = madeStep1, whichStep = whichStep1, checkDAGs = checkDAGs, TempState = TempState))
}


