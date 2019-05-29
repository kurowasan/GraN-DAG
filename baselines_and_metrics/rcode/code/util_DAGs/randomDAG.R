randomDAG <- function(p,probConnect,causalOrder = sample(p,p,replace=FALSE))
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
#    
# All rights reserved.  See the file COPYING for license terms. 
#
# simulates a directed acyclic graph (DAG) and returns its adjacency matrix
# 
# INPUT:
#   p           number of nodes 
#   probConnect the probability that an edge i -> j is added to the DAG
#   causalOrder starting with sink node (also called topological order)
#   
# OUTPUT:
#   DAG         Adjacency matrix of a directed acyclic graph (DAG)    
{
    #DAG <- as(diag(rep(0,p)),"sparseMatrix")
    DAG <- Matrix(nrow=p,ncol=p,0,sparse=TRUE)
    for(i in 1:(p-2))
    {
        node <- causalOrder[i]
        possibleParents <- causalOrder[(i+1):p]
        numberParents <- rbinom(n=1,size=(p-i),prob=probConnect)
        Parents <- sample(x = possibleParents, size = numberParents, replace = FALSE)
        DAG[Parents,node] <- rep(1,numberParents)
    }
    # Sample does not work properly when choosing from sets with one element. We thus consider the last case separately.  
    node <- causalOrder[p-1]
    ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
    DAG[causalOrder[p],node] <- ParentYesNo

    return(DAG)
}

scaleFreeDAG <- function(p,e,causalOrder = sample(p,p,replace=FALSE))
# Adapted from Jonas Peters' randomDAG function
# Use the Barabasi-Albert model
{
  DAG <- Matrix(nrow=p,ncol=p,0,sparse=TRUE)
  e <- e/p  # expected degree
  numberParents <- round(e, digits=0)
  
  # Sample does not work properly when choosing from sets with one element. We thus consider the last case separately.  
  origin <- causalOrder[1]
  node <- causalOrder[2]
  possibleParents <- causalOrder[1:2]
  DAG[origin,node] <- 1
  
  for(i in 3:p)
  {
    node <- causalOrder[i]
    Parents <- sample(x = possibleParents, size = numberParents, replace = TRUE) # TODO: make sure replace = TRUE
    DAG[Parents,node] <- rep(1,numberParents)
    
    possibleParents <- c(possibleParents, node)
    possibleParents <- c(possibleParents, Parents)
  }

  # returned the transposed of DAG, such that i->j
  return(t(DAG))  
}