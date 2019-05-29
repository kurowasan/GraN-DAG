iterAddNodes <- function(X, scoreName, pars, penFactor, output = FALSE)
# Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    computeNewState <- get(paste("computeNewState", scoreName, sep = ""))
    p <- dim(X)[2]
    initializeState <- get(paste("initialize",scoreName,sep=""))
    Vars <- sample(p,1)
    Perm <- 1
    SH <- abs(pars$SigmaHat) 
    diag(SH) <- rep(0,p)
    SH[Vars,] <- rep(0,p)
    
    for(i in 2:p)
    {
        if(1==1) #randomly select new nodes 
        {
            if(i<p)
            {
                new <- sample((1:p)[-Vars],1)
            }
            else
            {
                new <- (1:p)[-Vars]    
            }
        }
        else #select the node that has highest correlation with one of the exiting ones.
        {
            new <- which.max(SH[,Vars])%%p
            if(new ==0)
            {
                new <- p
            }
            SH[new,] <- rep(0,p)
        }
        Vars <- c(Vars,new)
        newPerm <- c(Perm,i)
        Xnew <- X[,Vars]
        
        initialAdj <- upper.tri(matrix(-1, i, i), diag = 0)
        if(scoreName == "SEMSEV")
        {
            pars <- computeStatisticsFromDataSEMSEV(Xnew,pars)
        }
        initialAdj[newPerm,newPerm] <- initialAdj
        newState <- initializeState(Xnew, initialAdj, pars)
        newState$CausOrder <- newPerm
        
        bestScore <- newState$Score
        bestPerm <- newPerm
        bestState <- newState
        
        for(j in i:2)
        {
            newState <- computeNewState(newState,c(j,j-1),Xnew,pars,penFactor,perm = 1)
            newPerm[c(j-1,j)] <- newPerm[c(j,j-1)] 
            if (newState$Score < bestScore)
            {
                bestScore <- newState$Score
                bestPerm <- newPerm
                bestState <- newState
            }        
        }
        Perm <-bestPerm
        if(output)
        {
            show(Vars[Perm])
        }
    }
    if(output)
    {
        cat("with score ", bestScore, " \n")
    }
    return(Vars[Perm])
}
