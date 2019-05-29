computeMinAdjSet <- function(AdjMat,i,j, pM)
    # Copyright (c) 2014 - 2014  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    #computes the minimal adjustment set from i to j. if j->i, we return just j
{
    p <- dim(AdjMat)[1]
    if(AdjMat[j,i] == TRUE)
    {
        return(j)
    }
    
    paSet <- which(AdjMat[,i]==TRUE) #choose the parents as current minimal set
    if(length(paSet)==0)
    {
        return(c())
    }
    
    #check whether emptyset is valid
    AdjMat2 <- AdjMat 
    AdjMat2[i,] <- rep(FALSE,p)
    if(dSepAdj(AdjMat2, i, j, c()))
    {
        return(c())
    }
    
    for(ii in 1:(length(paSet)-1))
    {
        allCombinations <- combinations(p,ii,1:p)
        numComb <- dim(allCombinations)[1]
        for(jj in 1:numComb)
        {
            currentSet <- allCombinations[jj,]
            
            if(isValidAdjSet(i,j,currentSet,AdjMat,pM))
            {
                return(currentSet)
            }
        }
    }
    
    return(paSet)
}