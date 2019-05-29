computePathMatrix <- function(G, spars=FALSE)
    # Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    #this function takes an adjacency matrix G from a DAG and computes a path matrix for which 
    # entry(i,j) being one means that there is a directed path from i to j
    # the diagonal will also be one
    p <- dim(G)[2]
    
    if((p > 3000) && (spars == FALSE))
    {
        warning("Maybe you should use the sparse version by using spars=TRUE to increase speed")
    }
    
    if(spars)
    {
        library(Matrix)
        G <- Matrix(G)
        PathMatrix <- Diagonal(p) + G
        
        k <- ceiling(log(p)/log(2))
        for(i in 1:k)
        {
            PathMatrix <- PathMatrix %*% PathMatrix            
        }
        PathMatrix <- PathMatrix > 0
        
    } else
    {
        # old and slow
        #PathMatrix <- G
        #    Gtmp <- diag(1,p)
        #    for(i in 1:p)
        #    {
        #        Gtmp <- Gtmp %*% G
        #        PathMatrix <- PathMatrix + Gtmp > 0
        #    }
        #    diag(PathMatrix) <- rep(1,p)    
        
        # new and faster
        PathMatrix <- diag(1,p) + G
        # PathMatrix <- as(diag(1,p) + G, "sparseMatrix")
        #   sparseMatrix does not seem to lead to a big improvement in speed (in fact, it seems slightly slower)
        
        k <- ceiling(log(p)/log(2))
        for(i in 1:k)
        {
            PathMatrix <- PathMatrix %*% PathMatrix            
        }
        
        PathMatrix <- PathMatrix > 0
    }
    return(PathMatrix)
}


