randomB <- function(G,lB = 0.1,uB = 0.9,twoIntervals = 1)
# if twoIntervals == TRUE, lB and uB should be positive
# Copyright (c) 2012-2012  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.
{
    numCoeff <- sum(G)
    B <- t(G)
    if(numCoeff ==1)
    {
        coeffs <- sample(c(-1,1),size=numCoeff,0.5)^(twoIntervals) * runif(1,lB,uB)
    }
    else
    {
        coeffs <- diag(sample(c(-1,1),size=numCoeff,0.5)^(twoIntervals)) %*% runif(numCoeff,lB,uB)
    }
    B[B==1] <- coeffs
    return(B)
}
