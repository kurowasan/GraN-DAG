connectedCompWrap <- function(G)
    # Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    conn.comp <- connectedComp(as(G,"graphNEL"))
    numConnComp <- length(conn.comp)
    for (ll in 1:numConnComp)
    {
        conn.comp[[ll]] <- as.numeric(conn.comp[[ll]])
    }
    return(conn.comp)
}