# Copyright (c) 2013 Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

set2n <- function(vec)
{
    if(length(vec)==0)
    {
        return(1)
    }
    else
    {
        return(sum(2^vec))
    }
}

allPartialCor <- function(corMat, output = FALSE)
    #requires library(gregmisc)
{
    p <- dim(corMat)[2]
    pCor <- array(rep(NA,p*p*2^p),dim = c(p,p,2^(p+1)))
    pCor[,,1] <- corMat
    for(i in 1:(p-2))
    {
        varsVec <- combinations(p, 2)
        for(numVars in 1:(dim(varsVec)[1]))
        {
            vars <- varsVec[numVars,]
            condSetVec <- combinations(p-2, i, (1:p)[-vars])
            for(numCondSet in 1:(dim(condSetVec)[1]))
            {
                condSet <- condSetVec[numCondSet,]    
                pCor[vars[1],vars[2],set2n(condSet)] <- (pCor[vars[1],vars[2],set2n(condSet[-1])] - pCor[vars[1],condSet[1],set2n(condSet[-1])] * pCor[condSet[1],vars[2],set2n(condSet[-1])]) / (sqrt(1-pCor[vars[1],condSet[1],set2n(condSet[-1])]^2) * sqrt(1-pCor[condSet[1],vars[2],set2n(condSet[-1])]^2))
                pCor[vars[2],vars[1],set2n(condSet)] <- pCor[vars[1],vars[2],set2n(condSet)]
                if(output)
                {
                    cat("partial correlation of ", vars[1], " and ", vars[2], " given ", condSet, " is ",pCor[vars[1],vars[2],set2n(condSet)], "\n")
                }
            }
        }
    }
    return(pCor)
}
