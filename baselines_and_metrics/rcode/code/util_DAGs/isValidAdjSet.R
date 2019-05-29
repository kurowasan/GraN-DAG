isValidAdjSet <- function(i,j,candSet,AdjMat,pM)
    # Copyright (c) 2014 - 2014  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    # checks whether a set is a valid adjustment set from i to j.
{
    PathMatrix2 <- computePathMatrix2(AdjMat,candSet,pM)
    
    p <- dim(AdjMat)[1]
    chiOfi <- which(AdjMat[i,]==TRUE)
    ancOfj <- which(pM[,j]==TRUE)
    
    chiCausPath <- intersect(chiOfi,ancOfj)
    if(length(chiCausPath)>0)
    {
        if(sum(pM[chiCausPath,candSet])>0)
        {
            return(FALSE)
        }
    }
    
    #check whether all non-causal paths are blocked
    checkAlldSep <- dSepAdji(AdjMat,i,candSet,pM,PathMatrix2)
    reachableWOutCausalPath <- checkAlldSep$reachableOnNonCausalPath
    if(reachableWOutCausalPath[j]==1)
    {
        return(FALSE)
    } else
    {
        return(TRUE)
    }
}