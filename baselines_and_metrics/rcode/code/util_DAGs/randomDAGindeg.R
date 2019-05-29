randomDAGindeg <- function(p,min_indeg,max_indeg, fillitup = FALSE)
    # Copyright (c) 2010 - 2013  	Jonas Peters [peters@stat.math.ethz.ch]
    #				                Jan Ernest [ernest@stat.math.ethz.ch] 
    # All rights reserved.  See the file COPYING for license terms. 
    
    # Input: p         - Number of vertices
    #        min_indeg - Minimum indegree of each vertex (0 or 1)
    #	     max_indeg - Maximum indegree of each vertex
    #        fillitup  - if set to FALSE, there are min_indeg roots. 
    #                    if set to TRUE, the last min_indeg nodes in the causal order are connected.
    
    # Output: DAG - random DAG satisfying min. and max. indegree constraints
{
    DAG <- diag(rep(0,p))
    causalOrder <- sample(p,p,replace=FALSE)  # starting with sink node
    
    for(i in 1:(p-(min_indeg+1)))
    {
        node <- causalOrder[i]
        possibleParents <- causalOrder[(i+1):p]
        if(max_indeg == min_indeg){
            numberParents <- min(min_indeg, p-i)
        } else {
            numberParents <- sample(min_indeg:min(max_indeg,p-i), size = 1)
        }	
        Parents <- sample(x = possibleParents, size = numberParents, replace = FALSE)
        DAG[Parents,node] <- rep(1,numberParents)
    }
    
    # As number of parents <= min_indeg, all remaining nodes are connected according to the causal order.  
    if(min_indeg > 0) 
    {
        if(fillitup)
        {
            for(j in 1:min_indeg)
            {
                node <- causalOrder[p-j]
                Parents <- causalOrder[(p-j+1):p]
                DAG[Parents,node] <- 1
            }
        } else
        {
            Parents <- causalOrder[(p-min_indeg+1):p]
            DAG[Parents, causalOrder[p-min_indeg]] <- 1
        }
    }
    return(DAG)
}

