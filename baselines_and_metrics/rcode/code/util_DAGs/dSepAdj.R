dSepAdj <- function(AdjMat, firstSet, secondSet, condSet)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    p <- dim(AdjMat)[2]
    A <- addDimNames(as.matrix(AdjMat))
    
    fir <- paste("x", firstSet, sep = "")
    sec <- paste("x", secondSet, sep = "")
    cond <- paste("x", condSet, sep = "")
    if(length(condSet) == 0)
    {
        return(msep(A,fir,sec,NULL))
    }
    else
    {
        return(msep(A,fir,sec,cond))
    }
        
}
