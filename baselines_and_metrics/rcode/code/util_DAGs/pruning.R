# Copyright (c) 2013-2013  Jonas Peters  [peters@stat.math.ethz.ch]
#                          Jan Ernest    [ernest@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
# 
# 
# Performs a pruning of the full DAG using a variable selection technique specified in pruneMethod.
# 
# INPUT:
#   X           nxp matrix of training inputs (n data points, p dimensions)
#   G           adjacency matrix of the full DAG
#   output      shall output be printed to the console
#   pruneMethod specifies the variable selection method. Default is selGam which fits a GAM with a specified
#               number of basis functions (numBasisFcts) vs all parents and selects using a threshold p-value
#               (cutOffPVal).
#   pruneMethodPars
#               supplies pruneMethod with parameter settings
#
# OUTPUT:
#   finalG      the pruned adjacency matrix    

pruning <- function(X, G, output = FALSE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10)) 
{
    p <- dim(G)[1]
    finalG <- matrix(0,p,p)
    for(i in 1:p)
    {
        parents <- which(G[,i]==1)
        lenpa <- length(parents)

        if(output)
        {
            cat("pruning variable:", i, "\n")
            cat("considered parents:", parents, "\n")
        }
        
        if(lenpa>0)
        {
            Xtmp <- cbind(X[,parents],X[,i])
            selectedPar <- pruneMethod(Xtmp, k = lenpa + 1, pars = pruneMethodPars, output = output)
            finalParents <- parents[selectedPar]
            finalG[finalParents,i] <- 1
        }
    }
    
    return(finalG)
}


# Performs a ranking of the edges according to their significance.
# 
# INPUT:
#   X           nxp matrix of training inputs (n data points, p dimensions)
#   G           adjacency matrix of the full DAG
#   output      shall output be printed to the console
#   rankMethod  specifies the ranking method. Default is rankGam.
#
# OUTPUT:
#   result      a list containing
#       $rankedEdges    a ranking of the edges in the Graph
#       $rankedPVal     a ranking of the pValues

rankingEdges <- function(X, G, output = FALSE, rankMethod = rankGam) 
{
    p <- dim(G)[1]
    numEdges <- sum(G)
    pValG <- matrix(Inf,p,p)
    for(i in 1:p)
    {
        if(output)
        {
            cat("\r Ranking performed for variable ", i, ".. \r")
        }
        parents <- which(G[,i]==1)
        lenpa <- length(parents)      
        if(lenpa>0)
        {
            Xtmp <- cbind(X[,parents],X[,i])
            pValPar <- rankMethod(Xtmp, k = lenpa + 1)
            pValG[parents,i] <- pValPar
        }
    }
    if(output)
    {
        cat("\n")
    }
    sss <- sort(pValG, index.return = T)
    rankedEdges <- arrayInd(sss$ix,c(p,p))
    rankedEdges <- rankedEdges[1:sum(G),]
    rankedPVal <- sss$x
    rankedPVal <- rankedPVal[1:sum(G)]
    result <- list(rankedEdges=rankedEdges, rankedPVal=rankedPVal)
    return(result)
}
