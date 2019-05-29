computeAndSaveAllFullDags <- function(p)
    # Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{   
    library(gtools)
    allPerms <- permutations(p,p)
    allFullDags <- matrix(NA,factorial(p),p^2)
    for(j in 1:factorial(p))
    {
        tmp <- matrix(0,p,p)
        for(i in 1:(p-1))
        {
            tmp[allPerms[j,i],allPerms[j,(i+1):p]] <- 1
        }
        allFullDags[j,] <- as.vector(tmp)
        #show(matrix(allFullDags[j,],p,p))
    }
    
    
    save(allFullDags, file = paste("~/Desktop/Link to svn/jpeters_causality/code/util_DAGs/allFullDagsWith",p,"Nodes", sep = ""))
    show(dim(allFullDags))
    
}
