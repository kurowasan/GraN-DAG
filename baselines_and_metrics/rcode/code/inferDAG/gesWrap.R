gesWrap <- function(X, lambda)
    # Copyright (c) 2012-2012  Jonas Peters [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms.
{
    result <- list()
    #score <- new("gauss.l0pen.obs.score", X)
    score <- new("GaussL0penObsScore", X, lambda=lambda)
    G <- ges(score)
    result$Adj <- as(G$essgraph, "matrix")
    n <- dim(X)[1]
    p <- dim(X)[2]
    #result$Score <- -(G$essgraph$score$global.score(G$repr) - n*p/2*log(2*pi) - n*p/2 )
    # JP 16.9.2013: I am not sure why this gives sth different than the next line 
    result$Score <- computeScoreSEMGauss(X,as(G$repr,"matrix"))
    return(result)
}

