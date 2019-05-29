computeNewStateSEMSEV <- function(oldState,index,X,pars,penFactor,perm=0)
    # Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
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
        #change the edge between i and j
    {
        newState$Adj[i,j] <- (newState$Adj[i,j] + 1) %% 2
        
        if(newState$Adj[j,i] == 1)
        {
            newState$Adj[j,i] <- 0
            recomputeNodes <- c(i,j)
        }
        else
        {
            recomputeNodes <- c(j)
        }
    }
    else
    {
        newState$CausOrder <- oldState$CausOrder
    }
    if(perm==1)
    {
        # exchange i and j
        perr <- 1:p
        perr[c(newState$CausOrder[i],newState$CausOrder[j])] <- perr[c(newState$CausOrder[j],newState$CausOrder[i])]
        newState$Adj <- newState$Adj[perr,perr]
        newState$CausOrder[c(i,j)] <- newState$CausOrder[c(j,i)]
        recomputeNodes <- newState$CausOrder[i:j]
#         show("we have to recompute")
#         show(recomputeNodes)
    }
    if(perm ==2)
    {
        # put i to position j 
        perr <- 1:p 
        tmpp <- newState$CausOrder
        tmpp2 <- rep(0,abs(i-j))
        if(i>j)
        {
            tmpp2 <- tmpp[j:(i-1)]            
            tmpp[j] <- tmpp[i]
            tmpp[(j+1):i] <- tmpp2
        }
        else
        {
            tmpp2 <- tmpp[(i+1):j]
            tmpp[j] <- tmpp[i]
            tmpp[i:(j-1)] <- tmpp2
        }
        
        perr <- oldState$CausOrder[sort(tmpp,index.return = TRUE)$ix]
        newState$CausOrder <- tmpp
        newState$Adj <- newState$Adj[perr,perr]
        
#        tmpi <- newState$CausOrder[i]
#        newState$CausOrder[j] <- tmpi
#        newState$Adj <- upper.tri(matrix(-1, p, p), diag = 0)
#        newState$Adj <- newState$Adj[newState$CausOrder,newState$CausOrder]
        
        recomputeNodes <- newState$CausOrder[i:j]
    }
    
    for(node in recomputeNodes)
    {
        newState$B[node,] <- rep(0,p)
        
        oldParents <- which(oldState$Adj[,node]==1)
        parents <- which(newState$Adj[,node]==1)
        
        # first, substract the part of the old variance
        newState$sumResVar = newState$sumResVar - newState$eachResVar[node]
        # and substract the number of old parents
        newState$numPars <- newState$numPars - length(oldParents)
        
        if(length(parents) == 0)
        {
            # use MLE instead of unbiased estimator of the variance
            newState$eachResVar[node] <- (n-1)/n*var(X[,node])
        }
        else
        {
            mod <- lm(X[,node] ~ X[,parents])  
            # use MLE instead of unbiased estimator of the variance
            newState$eachResVar[node] <- (n-1)/n*var(mod$residuals)
            newState$B[node,parents] <- mod$coef[2:(length(parents) + 1)]
        }
        # now, add the part of the new variance
        newState$sumResVar <- newState$sumResVar + newState$eachResVar[node]
        # and add the number of parameters
        newState$numPars <- newState$numPars + length(parents) # or + sum(State1$B[node,]!=0)        
    }
    
    newState$Score <- computeScoreSEMSEV(newState, pars, penFactor)
    return(newState)
}
