oneStepGPS <- function(StateOld, scoreName,pars,X,k,penFactor,checkDAGs)
    # Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    p <- StateOld$p
    # collect all different neighbors
    index <- array(dim = c(3*p*(p-1)/2,3))
    indexCount <- 0
    for(k in 1:min(p-1,15))
    {
        for(i in sample(p))
        {
            if(i+k <= p)
            {
                # moving one number 
                j <- i+k
                indexCount <- indexCount + 1
                index[indexCount,] <- c(i,j,2)
                
                # Transposition
                j <- i+k
                indexCount <- indexCount + 1
                index[indexCount,] <- c(i,j,1)
                

            }
            if(i-k >= 1)
            {
                # moving one number 
                j <- i-k
                indexCount <- indexCount + 1
                index[indexCount,] <- c(i,j,2)
            }            
        }
    }
    lenIndex <- indexCount
    # show(index)
    # permute this list
    #    index <- index[sample(lenIndex, replace = FALSE),]
    
    madeStep1 <- FALSE
    indexCount <- 0
    tried <- 0
    State1 <- StateOld
    whichStep1 <- c(-1,-1,-1)
    
    while( ((madeStep1 == FALSE)|(tried < k)) & (indexCount < lenIndex) )
    {    
        indexCount <- indexCount + 1
        i <- index[indexCount,1]
        j <- index[indexCount,2]
        per <- index[indexCount,3]
        
        tried <- tried +1
        checkDAGs <- checkDAGs + 1
        computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
        newState <- computeNewState(StateOld,c(i,j),X,pars,penFactor,perm=per)
        if(newState$Score < State1$Score)
        {
            State1 <- newState
            madeStep1 <- TRUE
            whichStep1 <- c(i,j,per)
            if(1==0)
            {
                show(c(i,j,per))
            }
        }
        
    }
    return(list(State = State1, madeStep = madeStep1, checkDAGs = checkDAGs, whichStep = whichStep1))
}


