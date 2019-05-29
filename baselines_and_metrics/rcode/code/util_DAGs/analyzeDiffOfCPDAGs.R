analyzeDiffOfCPDAGs <- function(trueG, estG)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    #both shoud be CPDAGs
    result <- list()
    p <- dim(trueG)[2]
    
    
    # skeleton
    trueGSkel <- pmin(trueG + t(trueG),1)
    estGSkel <- pmin(estG + t(estG),1)
    TP <- sum(trueGSkel * estGSkel)/2
    FP <- sum((1-trueGSkel) * estGSkel)/2
    P <- sum(trueGSkel)/2
    N <- p*(p-1)/2 - P
    result$TPRskel <- TP/P
    result$FPRskel <- FP/N
    
    # MEC
    TP <- sum(trueG * estG)
    FP <- sum((1-trueG) * estG)
    P <- sum(trueG)
    N <- p*(p-1)/2 - P
    result$TPR <- TP/P
    result$FPR <- FP/N
    
    return(result)
}


