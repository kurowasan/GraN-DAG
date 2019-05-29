lingamWrap <- function(X, output = FALSE)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    res <- lingam(t(X))
    return(list(B = res$B, Adj = t(res$B != 0)))
}
