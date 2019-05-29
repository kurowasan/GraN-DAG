standardizeMat <- function(X)
# Copyright (c) 2012-2012  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.
# standardizes the data (i.e. the columns have mean zero afterwards)
{
    if(is.null(dim(X)))
    {
        stop("X is not a matrix")
    }
    for(i in 1:dim(X)[2])
    {
        X[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])        
    }
    return(X)
}
