simulateFromRandomDAG <- function(p,probConnect,simulate)
    # Copyright (c) 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    DAG <- diag(rep(0,p))
    causalOrder <- sample(p,p,replace=FALSE)  # starting with sink node
    for(i in 1:(p-2))
    {
        node <- causalOrder[i]
        possibleParents <- causalOrder[(i+1):p]
        numberParents <- rbinom(n=1,size=(p-i),prob=probConnect)
        Parents <- sample(x = possibleParents, size = numberParents, replace = FALSE)
        DAG[Parents,node] <- rep(1,numberParents)
    }
    # Sample does not work properly when choosing from sets with one element. We thus consider the last case separately.  
    node <- causalOrder[p-1]
    ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
    DAG[causalOrder[p],node] <- ParentYesNo
    
    return(DAG,causalOrder,)
}
