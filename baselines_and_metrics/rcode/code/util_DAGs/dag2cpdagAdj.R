dag2cpdagAdj <- function(Adj)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    if(sum(Adj) == 0)
    {
        return(Adj)
    }
    cO <- computeCausOrder(Adj)
    d <- as(Adj[cO,cO], "graphNEL")
    cpd <- dag2cpdag(d)
    res <- matrix(NA,dim(Adj)[1],dim(Adj)[1])
    res[cO,cO] <- as(cpd, "matrix")
    result <- res
    ################
    #THE CODE ABOVE IS WRONG BECAUSE OF A VERY WEIRD BEHAVIOUR IN PCALG!!! Grrrr...
    ################
    ##Adj <- cbind(c(0,0,0,0),c(1,0,0,0),c(1,1,0,0),c(1,1,1,0))
    ##as(dag2cpdag(as(Adj, "graphNEL")), "matrix")
    ##plot(dag2cpdag(as(Adj, "graphNEL")))
    ##as(dag2cpdag(as(t(Adj), "graphNEL")), "matrix")
    ##plot(dag2cpdag(as(t(Adj), "graphNEL")))
    
    
    
    return(result)
}
