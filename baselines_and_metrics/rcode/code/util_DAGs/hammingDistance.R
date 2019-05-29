hammingDistance <- function(G1,G2, allMistakesOne = TRUE)
    # hammingDistance(G1,G2)
    #
    # Computes Hamming Distance between DAGs G1 and G2 with SHD(->,<-) = 1 if allMistakesOne == TRUE
    #
    # INPUT:  G1, G2     adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
    #         
    # OUTPUT: hammingDis Hamming Distance between G1 and G2
    #
    # Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms.
{
    if(allMistakesOne)
    {
        Gtmp <- (G1+G2)%%2
        Gtmp <- Gtmp + t(Gtmp)
        nrReversals <- sum(Gtmp == 2)/2
        nrInclDel <- sum(Gtmp == 1)/2
        hammingDis <- nrReversals + nrInclDel
    } else
    {        
        hammingDis <- sum(abs(G1 - G2))
        # correction: dist(-,.) = 1, not 2
        hammingDis <- hammingDis - 0.5*sum(G1 * t(G1) * (1-G2) * t(1-G2) + G2 * t(G2) * (1-G1) * t(1-G1))
    }    
    return(hammingDis)
}

computeFDR <- function(gt.DAG, exp.DAG)
    # fdr(G1,G2)
    #
    # Computes False Discovery Rate between estimated graph exp.DAG and ground truth DAG gt.DAG.
    # Supports only DAGs! not CPDAGs.
    #
    # INPUT:  gt.DAG, exp.DAG     ground truth DAG and estimated DAG
    #
    # OUTPUT: False discovery rate
    #
    # Written by authors of "Gradient-Based Neural DAG Learning"
{
    numerator <- sum((2 * exp.DAG + gt.DAG) == 2)
    nrEstimatedEdges <- sum(exp.DAG)
    return(numerator / nrEstimatedEdges)
}
