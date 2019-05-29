DiffOfDAGToPerm <- function(trueG, estG)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    #input: one DAG, one fully connected DAG!!
    
    result <- list()
    p <- dim(trueG)[2]
    
    dis <- sum(estG - trueG == -1)
    return(dis)
}
    
    
    
