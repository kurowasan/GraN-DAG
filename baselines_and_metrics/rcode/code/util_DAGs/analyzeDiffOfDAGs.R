analyzeDiffOfDAGs <- function(trueG, estG)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    #input: two DAGs!!
    
    result <- list()
    p <- dim(trueG)[2]
    
    # skeleton
    trueGSkel <- pmin(trueG + t(trueG),1) #pmin only for CPDAGs
    estGSkel <- pmin(estG + t(estG),1)
    TP <- sum(trueGSkel * estGSkel)/2
    FP <- sum((1-trueGSkel) * estGSkel)/2
    P <- sum(trueGSkel)/2
    N <- p*(p-1)/2 - P
    result$TPRskel <- TP/P
    result$FPRskel <- FP/N
    result$PRECskel <- TP/(TP+FP)
    
    # DAG
    TP <- sum(trueG * estG)
    FP <- sum((1-trueG) * estG)
    P <- sum(trueG)
    N <- p*(p-1) - P     # bidirected edge = 2 edges
    result$TPR <- TP/P
    result$FPR <- FP/N
    result$PREC <- TP/(TP+FP)
    
    return(result)
}
    
    
    
