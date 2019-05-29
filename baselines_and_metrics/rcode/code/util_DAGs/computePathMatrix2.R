computePathMatrix2 <- function(G,condSet,PathMatrix1)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    # The only difference to the function computePathMatrix is that this function changes
    # the graph by removing all edges that leave condSet.
    # If condSet is empty, it just returns PathMatrix1.
    
    p <- dim(G)[2]
    if(length(condSet) > 0)
    {
        G[condSet,]<-matrix(0,length(condSet),p)
        
        
        # old and slow
        #PathMatrix2 <- G
        #Gtmp <- diag(1,p)
        #    
        #for(i in 1:p)
        #{
        #    Gtmp <- Gtmp %*% G
        #    PathMatrix2 <- PathMatrix2 + Gtmp > 0
        #}
        #diag(PathMatrix2) <- rep(1,p)    
        
        # new and faster
        PathMatrix2 <- diag(1,p) + G
        #PathMatrix2 <- as(diag(1,p) + G, "sparseMatrix")
        k <- ceiling(log(p)/log(2))
        for(i in 1:k)
        {
            PathMatrix2 <- PathMatrix2 %*% PathMatrix2            
        }
        PathMatrix2 <- PathMatrix2 > 0        
    } else
    {
        PathMatrix2 <- PathMatrix1
    }
    return(PathMatrix2)
}
    
