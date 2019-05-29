# Copyright (c) 2014 - 2014  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

allDagsJonas <- function (adj,row.names)
{
    # Input: adj. mat of a DAG with row.names, containing the undirected component that 
    #        should be extended
    # !!!! the function can probably be faster if we use partial orderings 
    
    a <- adj[row.names,row.names]
    
    if(any((a + t(a))==1))
    {
        #warning("The matrix is not entirely undirected.")
        return(-1)
    }
    return(allDagsIntern(adj, a, row.names,NULL))
}

allDagsIntern <- function (gm, a, row.names, tmp) 
    # Input: adj. mat gm, submatrix a of an UNDIRECTED component and tmp that is set to NULL 
    # (due to the recursive nature of the function),     
{
    if(any((a + t(a))==1))
    {
        stop('The matrix is not entirely undirected. This should not happen!')
    }
    
    if (sum(a) == 0) 
    {
        tmp2 <- rbind(tmp, c(gm))
        if (all(!duplicated(tmp2))) # if tmp2 contains no element twice
            tmp <- tmp2
    }
    else 
    {
        # all nodes can be sinks, but we consider only those who have neighbors.
        sinks <- which(colSums(a) > 0)
        for (x in sinks) 
        {
            gm2 <- gm
            a2 <- NULL
            row.names2 <- NULL
            
            Adj <- (a == 1)
            Adjx <- Adj[x,]
            if(any(Adjx))
            {
                un <- which(Adjx)
                pp <- length(un)
                Adj2 <- matrix(Adj[un,un],pp,pp)
                diag(Adj2) <- rep(TRUE,pp)
            } else # x does not have any neighbors
            {
                Adj2 <- TRUE
            }
            # Are all (undirected) neighbors of x connected? (O/wise there will be a v-structure if 
            # x becomes a sink node)
            if(all(Adj2)) #if not, don't do anything
            {
                if(any(Adjx))
                {
                    un <- row.names[which(Adjx)]
                    pp <- length(un)
                    gm2[un,row.names[x]] <- rep(1,pp)
                    gm2[row.names[x],un] <- rep(0,pp)
                }
                a2 <- a[-x, -x]
                row.names2 <- row.names[-x]
                tmp <- allDagsIntern(gm2, a2, row.names2, tmp)
            }
        }
    }
    return(tmp)
}
