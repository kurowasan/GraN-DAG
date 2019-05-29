dSepAdji <- function(AdjMat, i, condSet, PathMatrix = computePathMatrix(AdjMat), PathMatrix2 = matrix(NA,p,p))
    # Copyright (c) 2013 - 2014  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    # This function looks for all j, such that i is d-sep to j given the condSet
    # The PathMatrix contains ancestor-relations.
    #
    # REACHABLE: 
    #      means that there is a path from i that is not blocked by condSet 
    # REACHABLEONNONCAUSALPATH = REACHABLEONNONDIRECTEDPATH (both notations exist):
    #      means that there is a path from i that is not blocked by condSet and that is NOT directed
    #      
    # PathMatrix2 contains ancestor-relations when removing arrows condSet->
    # reachableOnNonCausalPath indicates whether in AdjMat there is a nonDirected path from i to j that is not blocked by condSet. 
    # More details: For all combinations of nodes and edge direction (incoming or outgoing), we check, to which 
    # neighbouring node-orientation pair an open path can be propagated
    # further. The main trick is to construct a $2p \times 2p$ matrix. The value $p+j$, for example, stands for
    # node $j$ with an outgoing edge. If the entry $(j_1,p+j_2)$ is one, it means that if we can reach
    # $j_1$ by $\rightarrow j_1$ from node $i$, then we can also reach $j_2$ with a path ending $\leftarrow j_2$. 
{
    if(is.na(sum(PathMatrix2)))
    {
        PathMatrix2 <- computePathMatrix2(AdjMat,condSet,PathMatrix)    
    }
    if(length(condSet)==0)
    {
        AncOfCondSet <- c()
    }
    if(length(condSet)==1)
    {
        AncOfCondSet <- which(PathMatrix[,condSet] > 0)
    }
    if(length(condSet) > 1)
    {
        AncOfCondSet <- which(rowSums(PathMatrix[,condSet])>0)
    }
    p <- dim(AdjMat)[2]
    reachabilityMatrix <- matrix(0,2*p,2*p)
    reachableOnNonCausalPathLater <- matrix(0,2,2)
    
    reachableNodes <- rep(0,2*p) # the first p entries contain the nodes that are reachable using an incoming edge (e.g. children of i), the last p entries the ones reachable using an outgoing edge (e.g. parents of i)
    reachableOnNonCausalPath <- rep(0,2*p) # same idea as above
    alreadyChecked <- rep(0,p)
    k <- 2
    toCheck <- c(0,0)
    
    reachableCh <- which(AdjMat[i,] == 1)
    if(length(reachableCh) > 0)
    {
        toCheck <- c(toCheck,reachableCh)
        reachableNodes[reachableCh] <- rep(1,length(reachableCh)) 
        AdjMat[i,reachableCh] <- rep(0,length(reachableCh))
    }
    
    reachablePa <- which(AdjMat[,i] == 1)
    if(length(reachablePa) > 0)
    {
        toCheck <- c(toCheck,reachablePa)
        reachableNodes[reachablePa+p] <- rep(1,length(reachablePa)) 
        reachableOnNonCausalPath[reachablePa + p] <- rep(1,length(reachablePa))
        AdjMat[reachablePa,i] <- rep(0,length(reachablePa))
    }
    
    
    while(k < length(toCheck))
    {
        k <- k + 1
        a1 <- toCheck[k]
        if(alreadyChecked[a1] == 0)
        {
            #  -> currentNode => a2 == 1
            currentNode <- a1
            alreadyChecked[a1] <- 1
            
            
            ######
            #PARENTS OF CURRENTNODE
            ######
            Pa <- which(AdjMat[,currentNode] == 1)
            
            # IF   one of the Pa of currentNode is reachable and is not included in condSet, 
            # THEN currentNode is, too
            Pa1 <- setdiff(Pa,condSet)
            reachabilityMatrix[Pa1,currentNode] <- rep(1,length(Pa1))
            reachabilityMatrix[Pa1 + p,currentNode] <- rep(1,length(Pa1))
            
            # IF   currentNode is reachable with -> cN and cN is in AncOfCondSet, 
            # THEN parents are, too
            if(sum(AncOfCondSet == currentNode)>0)
            {
                reachabilityMatrix[currentNode,Pa + p] <- rep(1,length(Pa))  #not necessary?
                if(PathMatrix2[i,currentNode] > 0)
                {
                    reachableOnNonCausalPathLater <- rbind(reachableOnNonCausalPathLater,cbind(rep(currentNode,length(Pa)),Pa))
                }
                newtoCheck <- Pa
                # remove nodes that have already been checked
                newtoCheck <- newtoCheck[which(alreadyChecked[newtoCheck] == 0)]
                toCheck <- c(toCheck,newtoCheck)
            }
            
            # IF   currentNode is reachable with <- cN and cN is in not in condSet, 
            # THEN parents are, too
            if(sum(condSet == currentNode) == 0)
            {
                reachabilityMatrix[currentNode + p,Pa + p] <- rep(1,length(Pa))    
                newtoCheck <- Pa
                # remove nodes that have already been checked
                newtoCheck <- newtoCheck[which(alreadyChecked[newtoCheck] == 0)]
                toCheck <- c(toCheck,newtoCheck)
            }
            
            
            ######
            #CHILDREN OF CURRENTNODE
            ######
            Ch <- which(AdjMat[currentNode,] == 1)
            
            # IF    Ch of currentNode is reachable on a path with <- Ch and Ch is not in CondSet
            # THEN  currentNode is reachable, too
            Ch1 <- setdiff(Ch, condSet)
            reachabilityMatrix[Ch1 + p, currentNode + p] <- rep(1,length(Ch1))    
            
            # IF Ch of currentNode is reachable on a path with -> Ch and Ch is in AncOfCondSet
            # THEN  currentNode is reachable, too
            Ch2 <- intersect(Ch,AncOfCondSet)
            reachabilityMatrix[Ch2, currentNode + p] <- rep(1,length(Ch2))   #not necessary? 
            Ch2b <- intersect(Ch2,which(PathMatrix2[i,] > 0))
            reachableOnNonCausalPathLater <- rbind(reachableOnNonCausalPathLater,cbind(Ch2b,rep(currentNode,length(Ch2b))))
            
            # IF   currentNode is reachable and cN is in not in condSet, 
            # THEN Ch are, too
            if(sum(condSet == currentNode) == 0)
            {
                reachabilityMatrix[currentNode,Ch] <- rep(1,length(Ch))    
                reachabilityMatrix[currentNode+p,Ch] <- rep(1,length(Ch))  
                newtoCheck <- Ch
                # remove nodes that have already been checked
                newtoCheck <- newtoCheck[which(alreadyChecked[newtoCheck] == 0)]
                toCheck <- c(toCheck,newtoCheck)
            }
        }
    }
    
    #reachabilityMatrix <- as(reachabilityMatrix, "sparseMatrix")
    reachabilityMatrix <- computePathMatrix(reachabilityMatrix)
    reachabilityMatrix <- as(reachabilityMatrix, "matrix")
    #that is expensive!!!
    
    
    ttt2 <- which(reachableNodes == 1)
    if(length(ttt2)==1)
    {
        tt2 <- which(reachabilityMatrix[ttt2,] > 0)
    } else
    {
        tt2<-which(colSums(reachabilityMatrix[ttt2,])>0)        
    }
    reachableNodes[tt2]<-rep(1,length(tt2))
    
    
    
    #first activation step (all nodes that are reachable from the parents of i are reachable on a nonDirected = nonCausal path)
    ttt <- which(reachableOnNonCausalPath == 1)
    if(length(ttt)==1)
    {
        tt <- which(reachabilityMatrix[ttt,] > 0)
    } else
    {
        tt<-which(colSums(reachabilityMatrix[ttt,])>0)        
    }
    reachableOnNonCausalPath[tt]<-rep(1,length(tt))
    
    #second activation step
    # there is one type of nodes that are reachable on a nondirected path that we have missed so far:
    # nodes that are reachable from parents of k where k is (a descendant of i and reachable from i) 
    # the matrix reachableOnNonCausalPathLater already contains reachable parents of reachable nodes. They are therefore reachable on a nonCausal path.
    if(dim(reachableOnNonCausalPathLater)[1] > 2)
    {
        for(kk in 3:(dim(reachableOnNonCausalPathLater)[1]))
        {
            ReachableThrough <- reachableOnNonCausalPathLater[kk,1]        
            newReachable <- reachableOnNonCausalPathLater[kk,2]
            reachableOnNonCausalPath[newReachable+p] <- 1
            
            # cancel the connection from the parent of the reachable node (newReachable) to the reachable node (ReachableThrough); otherwise  ReachableThrough seems to be reachable on a nonCausal path
            reachabilityMatrix[newReachable,ReachableThrough] <- 0
            reachabilityMatrix[newReachable,ReachableThrough+p] <- 0
            reachabilityMatrix[newReachable+p,ReachableThrough] <- 0
            reachabilityMatrix[newReachable+p,ReachableThrough+p] <- 0
        }
        ttt <- which(reachableOnNonCausalPath == 1) # we only have to check the newReachables but we check all reachableOnNonCausalPath again; avoiding this may save some time
        if(length(ttt)==1)
        {
            tt <- which(reachabilityMatrix[ttt,] > 0)
        } else
        {
            tt<-which(colSums(reachabilityMatrix[ttt,])>0)        
        }
        reachableOnNonCausalPath[tt]<-rep(1,length(tt))
    }
    
    # show(reachableOnNonCausalPath)
    # show(toCheck)
    result <- list()
    result$reachableJ <- rowSums(cbind(reachableNodes[1:p],reachableNodes[(p+1):(2*p)])) > 0
    result$reachableOnNonCausalPath <- rowSums(cbind(reachableOnNonCausalPath[1:p],reachableOnNonCausalPath[(p+1):(2*p)])) > 0
    return(result)
}
