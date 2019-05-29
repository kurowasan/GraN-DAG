computeSigma <- function(B,noisevar)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    p <- dim(B)[2]
    Id <- diag(rep(1,p))
    hilfsm <- solve(Id-B) 
    SigmaX <- hilfsm %*% diag(noisevar) %*% t(hilfsm)
    return(SigmaX)
}
