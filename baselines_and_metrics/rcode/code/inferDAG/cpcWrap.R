cpcWrap <- function(X, alphaU, mmax)
    # Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms.
{
    suffStat <- list(C = cor(X), n = nrow(X))
    pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(X), alpha = alphaU, m.max = mmax, conservative = TRUE, u2pd = "relaxed")
    result <- as(pc.fit@graph, "matrix")
    return(result)
}
