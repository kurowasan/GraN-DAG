analyzeCPDAG <- function(A)
# Copyright (c) 2012-2012  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.
{
    li <- list()
    li$numEdges <- sum((A + t(A))>0)/2     
    li$numUndirectedEdges <- sum(A) - li$numEdges    
    li$numDirectedEdges <- li$numEdges -  li$numUndirectedEdges    
    return(li)
}
