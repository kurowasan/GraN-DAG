addDimNames <- function(AdjMat)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    p <- dim(AdjMat)[2]
    formul <- paste("x1 ~ x2", sep = "")
    for(i in 3:p)
    {
        formul <- paste(formul, " + x",i, sep = "")
    }
    A <- DAG(formula(formul))
    A[,] <- AdjMat
    return(A)
}
