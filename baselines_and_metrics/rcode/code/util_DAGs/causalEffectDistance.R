causalEffectDistance <- function(G1,G2,SigmaX)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    p <- dim(G1)[2]
    cE1 <- matrix(0,p,p)
    cE2 <- matrix(0,p,p)
    for(i in 1:p)
    {
        for(j in 1:p)
        {
            if(i!=j)
            {
                cE1[i,j] <- computeCausalEffect(i,j,SigmaX,G1)
                cE2[i,j] <- computeCausalEffect(i,j,SigmaX,G2)
            }
        }
    }
    return(sum(abs(cE1 - cE2)>10^(-8)))
}
