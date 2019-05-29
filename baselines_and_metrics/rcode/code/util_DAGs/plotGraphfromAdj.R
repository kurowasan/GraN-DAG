plotGraphfromAdj <- function(Adj, is_weighted = NULL)
    # Copyright (c) 2013 - 2013  Jan Ernest [ernest@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    
    # Input:  Adj         - Adjacency matrix of a graph
    #         is_weighted - TRUE if the plot should contain edge weights
    # Output: Plot of the corresponding (un-)weighted directed graph
    
{
    library(igraph)
    
    if(!is.null(is_weighted)) 
    { 
        gr <- graph.adjacency(Adj, mode = "directed", weighted = is_weighted, diag = FALSE)	
        E(gr)$label <- E(gr)$weight 
    } else 
    {
        Adj[Adj != 0] <- 1 
        gr <- graph.adjacency(Adj, mode = "directed", weighted = is_weighted, diag = FALSE)
#        gr <- graph.adjacency(Adj, mode = "undirected", weighted = is_weighted, diag = FALSE)
    }
    
    V(gr)$label <- V(gr) 
    
    plot(gr)
}
