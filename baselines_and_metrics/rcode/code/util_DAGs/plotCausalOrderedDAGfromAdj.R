plotCausalOrderedDAGfromAdj <- function(Adj, labels = 1:dim(Adj)[1], main=NULL)
    # Copyright (c) 2013 - 2013  Jan Ernest [ernest@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    
    # Input:  Adj         - Adjacency matrix of a graph
    #         labels      - an optional vector of labels of the nodes in the graph
    #         
    # Output: Plot of the corresponding (un-)weighted directed graph respecting the causal order. 
    
{
    library(Rgraphviz)
    G <- as(Adj, "graphNEL")
    
    z <- labels
    names(z) = nodes(G)
    nAttrs <- list()
    nAttrs$label <- z
    attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE))
    plot(G, nodeAttrs = nAttrs, attrs = attrs, main=main)    
}