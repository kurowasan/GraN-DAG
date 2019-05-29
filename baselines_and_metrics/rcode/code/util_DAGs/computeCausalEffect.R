computeCausalEffect <- function(i,j,SigmaX,AdjMat)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
#computes the causal Effect from i to j in the case of Gaussian variables
{
    if(AdjMat[j,i] == 1)
    {
        res <- 0
    } else
    {
        pai <- which(AdjMat[,i] == 1)
        regressors <- c(i,pai)
        Sigma22 <- SigmaX[regressors, regressors]
        Sigma12 <- SigmaX[j,regressors]
        coefs <- Sigma12 %*% solve(Sigma22)
        res <- coefs[1]
    }
    
    #pruning
#     if(abs(res) < 10^(-12))
#     {
#         res <- 0
#     }
    return(res)
}
