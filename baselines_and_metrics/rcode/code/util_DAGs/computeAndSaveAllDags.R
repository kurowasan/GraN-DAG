computeAndSaveAllDags <- function(p)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{   
    allExpansionsWoutDiag <- as.matrix(expand.grid( rep(list(0:1),p*(p-1)) ))
    numCandidates <- dim(allExpansionsWoutDiag)[1]
    allCandidates <- matrix(0,numCandidates,p*p)
    for(j in 1:(p-1))
    {
        # which stays zero? j+(j-1)*p
         allCandidates[,(j+(j-1)*p+1):((j+1)+j*p -1)] <- allExpansionsWoutDiag[,((j-1)*p+1):(j*p)]  
    }
    
    validDag <- rep(TRUE,numCandidates)
    cat("starting now \n")
    for(i in 1:numCandidates)
    {
        cat(i," out of ",numCandidates, "\r")
        G <- matrix(allCandidates[i,],p,p)
        if(containsCycle(G))
        {
            validDag[i] <- FALSE
        }
    }
    cat("\n DONE! \n")
    allDags <- allCandidates[which(validDag),]
    allDags <- (allDags == 1)
    save(allDags, file = paste("allDagsWith",p,"Nodes", sep = ""))
    show(dim(allDags))
    
}
