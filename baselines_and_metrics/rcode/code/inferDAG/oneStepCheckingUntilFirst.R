oneStepCheckingUntilFirst <- function(State1, scoreName,pars,X,penFactor,checkDAGs)
    # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    TempState <- State1
    TempState$Score <- Inf
    TempState$ScoreStat <- Inf
    p <- State1$p
    
    if(scoreName == "SEMSEV")
    {
        # start with largest residual variance.
        # sortNodes <- sort(StateOld$eachResVar, index.return=TRUE, decreasing=TRUE)$ix
        varTmp <- State1$eachResVar - min(State1$eachResVar)
        #varTmp <- abs(StateOld$eachResVar - mean(StateOld$eachResVar))
        varTmp <- varTmp + min(varTmp[varTmp>0])
        varTmp <- varTmp / sum(varTmp)
        sortNodes <- sample(x = 1:p, size = p, p = varTmp, replace = FALSE)
    } else
    {
        if(scoreName == "SEMIND")
        {
            sortNodes <- sample(x = 1:p, size = p, p = pmax(State1$SingleScores,0.000001), replace = FALSE)   
        } else
        {
            sortNodes <- 1:p
        }
    }
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
    # show(State1$SingleScores)
    # show(index)
    # permute this list
    # index <- index[sample(p*(p-1), replace = FALSE),]
    
    madeStep1 <- FALSE
    indexCount <- 0
    
    while( (madeStep1 == FALSE) & (indexCount < (p*(p-1))) )
    {    
        indexCount <- indexCount + 1
        i <- index[indexCount,1]
        j <- index[indexCount,2]
        
        candidateAdj <- State1$Adj
        candidateAdj[i,j] <- (candidateAdj[i,j] + 1) %% 2
        candidateAdj[j,i] <- 0
        
        if((!containsCycle(candidateAdj)) & ((i != State1$tabuInd[1]) | (j != State1$tabuInd[2])) )
        {
            checkDAGs <- checkDAGs + 1
            computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
            newState <- computeNewState(State1,c(i,j),X,pars,penFactor)
            
            if( (scoreName == "SEMIND") && (max(State1$SingleScores) == 100) )
            {
                if(newState$ScoreStat < State1$ScoreStat)
                {
                    show(State1$ScoreStat)
                    show(newState$ScoreStat)
                    State1 <- newState
                    State1$tabuInd <- c(i,j)
                    madeStep1 <- TRUE                    
                }
                if(newState$ScoreStat < TempState$ScoreStat)
                {
                    TempState <- newState
                    TempState$tabuInd <- c(i,j)
                }
            } else
            {
                if(newState$Score < State1$Score)
                {
                    #show(State1$Score)
                    #cat("Hallo1a\n ")
                    #show(newState$Score)
                    #cat("Hallo1b\n ")
                    State1 <- newState
                    State1$tabuInd <- c(i,j)
                    madeStep1 <- TRUE
                }
                if(newState$Score < TempState$Score)
                {
                    TempState <- newState
                    TempState$tabuInd <- c(i,j)
                }
            }
        }
    }
    TempState$Score <- State1$Score - 10^(-30)
    
    return(list(State = State1, madeStep = madeStep1, checkDAGs = checkDAGs, TempState = TempState))
}


