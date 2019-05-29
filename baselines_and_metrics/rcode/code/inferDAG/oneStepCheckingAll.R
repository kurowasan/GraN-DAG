oneStepCheckingAll <- function(State,scoreName,pars,X,penFactor,checkDAGs)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    p <- State$p
    bestTillNowState <- State
    madeStep1 <- FALSE
    for(i in 1:p)
    {
        for(j in 1:(p-1))
        {
            # do not add sth on the diagonal
            j <- j + (j>(i-1))
            index <- c(i,j)
                
            candidateAdj <- State$Adj
            candidateAdj[i,j] <- (candidateAdj[i,j] + 1) %% 2
            candidateAdj[j,i] <- 0
                    
            if(!containsCycle(candidateAdj))
            {
                checkDAGs<- checkDAGs + 1
                computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
                newState <- computeNewState(State,index,X,pars,penFactor)
                if(scoreName == "SEMINDDDDDDD")
                {
                    if(newState$DiffScore < bestTillNowState$DiffScore)
                    {
                        bestTillNowState <- newState
                        madeStep1 <- TRUE
                    }
                }
                else
                {
                    if(newState$Score < bestTillNowState$Score)
                    {
                        bestTillNowState <- newState
                        madeStep1 <- TRUE
                    }
                }
                
            }
        }        
    }
    State1 <- bestTillNowState
    State1$DiffScore <- 0
    return(list(State = State1, madeStep = madeStep1, checkDAGs = checkDAGs))
}


