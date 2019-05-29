containsCycle <- function(Adj)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    ifelse(sum(abs(eigen(Adj, symmetric = FALSE, only.values = TRUE)$values) > 10^(-10)) > 0, TRUE, FALSE)    
}
