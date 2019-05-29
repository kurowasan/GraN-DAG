structIntervDist <- function(trueGraph, estGraph, output = FALSE)
    # Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    # These are all initializations
    estGraph <- as(estGraph, "matrix") #to make the code better readable we sometimes write Gp instead of estGraph
    trueGraph <- as(trueGraph, "matrix") #to make the code better readable we sometimes write G instead of trueGraph
    p <- dim(trueGraph)[2]
    rownames(estGraph) <- 1:p
    colnames(estGraph) <- 1:p    
    incorrectInt <- matrix(0,p,p)
    correctInt <- matrix(0,p,p)
    minimumTotal <- 0
    maximumTotal <- 0
    timePathMatrix2 <- 0    
    timeAlldSep <- 0    
    timeExpGraph <- 0
    numChecks <- 0
    Gtmp <- diag(1,p)
    ptmtotal <- proc.time()
    
    # Compute the path matrix whose entry (i,j) is TRUE if there is a directed path
    # from i to j. The diagonal is TRUE, too.
    PathMatrix <- computePathMatrix(trueGraph)
    PathMatrix <- as(PathMatrix, "matrix")
    
    # We now compute the undirected graph and all its connected components.
    # The graph does not contain undirected components if it is a DAG.  
    Gp.undir <- estGraph * t(estGraph)
    gp.undir <- as(Gp.undir,"graphNEL")
    conn.comp <- connectedComp(gp.undir)
    numConnComp <- length(conn.comp) # == p if (but not only if) estGraph is a DAG
    GpIsEssentialGraph <- TRUE    
    for (ll in 1:numConnComp) 
    {
        conn.comp[[ll]] <- as.numeric(conn.comp[[ll]])
        if(length(conn.comp[[ll]]) > 1)
        {
            chordal <- is.chordal(igraph.from.graphNEL(as(Gp.undir[conn.comp[[ll]],conn.comp[[ll]]],"graphNEL")),fillin=TRUE)$chordal
            if(!chordal)
            {
                show("The estimated graph is not chordal, i.e. it is not a CPDAG! We thus consider local expansions of the graph (some combinations of which may lead to cycles).")
                GpIsEssentialGraph <- FALSE
            } 
            if(length(conn.comp[[ll]]) > 8)
            {
                show("The connected component is too large (>8 nodes) in order to be extended to all DAGs in a reasonable amount of time. We thus consider local expansions of the graph (some combinations of which may lead to cycles).")
                GpIsEssentialGraph <- FALSE
            }          
        }
    }
    
    
    for (ll in 1:numConnComp)
    {
        ptm <- proc.time()
        if(length(conn.comp[[ll]]) > 0)
        {
            if(GpIsEssentialGraph)
            {
                # expand the connected component into DAGs
                mmm <- allDagsJonas(estGraph,conn.comp[[ll]])
                if(sum(mmm == -1) == 1)
                {
                    GpIsEssentialGraph <- FALSE
                    mmm <- matrix(c(estGraph),1,p^2)
                }
                if(dim(mmm)[1] == 1)
                {
                    if(sum( abs( matrix(mmm,p,p) - estGraph ) ) > 0)
                    {
                        stop("WHAT????")
                    }
                }
                # each row in mmm contains one DAG.
                # the first p entries of estGraph are the parents of node 1
                # we change this to: the first p entries of estGraph are the children of node 1
                newInd <- seq(1,p^3,by=p)-rep(seq(0,(p-1)*(p^2-1),by=(p^2-1)),each=p)
                dimM <- dim(mmm)
                mmm <- matrix(mmm[,newInd],dimM)
                
                if(is.null(mmm))
                {
                    show("Something is wrong. Maybe the estimated graph is not a CPDAG? We expand the undirected components locally.")
                    GpIsEssentialGraph <- FALSE
                } else 
                {
                    incorrectSum <- rep(0,dim(mmm)[1])
                }
            } 
        }
        timeExpGraph <- timeExpGraph + proc.time() - ptm
        
        for(i in conn.comp[[ll]]) # this for loop is over nodes in the conn. comp.
        {   
            paG <- which(trueGraph[,i]==1) #parents of i in trueGraph             
            certainpaGp <- which((estGraph[,i]*(rep(1,p)-estGraph[i,])) == 1) # these nodes are parents of i in estGraph
            possiblepaGp <- which((estGraph[,i]*estGraph[i,]) == 1) # all nodes j s.t. (i,j) == (j,i) == 1 in estGraph 
            if(!GpIsEssentialGraph) #in this case go through all local combinations of parents. We do not care whether this is consistent with a graph structure.
            {
                maxcount <- 2^length(possiblepaGp) #nr of possible parent sets
                uniqueRows <- 1:maxcount
                mmm <- t(matrix(rep(t(estGraph)[1:length(estGraph)],maxcount),length(estGraph),maxcount)) #write them into mm (each row is a matrix)
                mmm[,i+(possiblepaGp-1)*p] <- as.matrix(expand.grid( rep( list(0:1), length(possiblepaGp) ) )) #cont.
                incorrectSum <- rep(0,maxcount) #for each i this will become the nr of j s.t. the intervention (i,j) is incorrect
            } else 
            {   
                if(dim(mmm)[1]>1)
                {
                    # each row in mmm contains a different DAG expansion of the l-th conn.comp. However, the parent sets of node i might be the same for many DAGs. 
                    # construct uniqueRows that indicate the rows with different sets of parents for node i.
                    allParentsOfI <- seq(i,(p-1)*p+i,by=p)
                    uniqueRows <- which(!duplicated(mmm[,allParentsOfI])) #of mmm
                    maxcount <- length(uniqueRows) #maxcount is the number of different parent sets for node i in the ll-th conn. comp.
                } else #this means we have a DAG
                {
                    maxcount <- 1
                    uniqueRows <- 1
                }
            }
            
            count <- 1
            while(count <= maxcount) # go through all expansions. In case of a CPDAG, we only go through uniqueRows which indicate the unique set of parents of i.
            {
                if(maxcount == 1) #maxcount is the number of (local) expansions we consider. if it is 1, we have only one set of possible parents
                {
                    paGp <- certainpaGp
                } else
                {
                    Gpnew <- t(matrix(mmm[uniqueRows[count],],p,p))
                    paGp <- which(Gpnew[,i] == 1) # the set of parents we have to check
                    if(output)
                    {
                        cat(i," has ", length(paGp), " parents in expansion nr. ", uniqueRows[count], " of Gp:")
                        show(paGp)
                    }
                }
                
                #the following computations are the same for all j (i is fixed)
                ptm <- proc.time()
                PathMatrix2 <- computePathMatrix2(trueGraph,paGp,PathMatrix)
                timePathMatrix2 <- timePathMatrix2 + proc.time() - ptm
                
                ptm <- proc.time()
                checkAlldSep <- dSepAdji(trueGraph,i,paGp,PathMatrix,PathMatrix2)
                numChecks <- numChecks + 1
                timeAlldSep <- timeAlldSep + proc.time() - ptm
                
                reachableWOutCausalPath <- checkAlldSep$reachableOnNonCausalPath
                
                
                for(j in 1:p)
                {
                    if(i != j) # test the intervention effect from i to j
                    {
                        # The order of the following checks and the flag finished are
                        # made such that as few tests are performed as possible. 
                        
                        finished <- FALSE
                        ijGNull <- FALSE
                        ijGpNull <- FALSE
                        
                        # ijGNull means that the causal effect from i to j is zero in G
                        # more precisely, p(x_j | do (x_i=a)) = p(x_j) 
                        if(PathMatrix[i, j] == 0)
                        {
                            ijGNull <- TRUE
                        }               
                        
                        # if j->i exists in Gp
                        if( (sum(paGp==j)==1) )
                        {
                            ijGpNull <- TRUE
                        }
                        
                        # if both are zero
                        if(ijGpNull & ijGNull)
                        {
                            finished <- TRUE
                            correctInt[i,j] <- 1
                        }
                        
                        # if Gp predicts zero but G says it is not
                        if(ijGpNull & !ijGNull)
                        {
                            incorrectInt[i,j] <- 1
                            incorrectSum[uniqueRows[count]] <- incorrectSum[uniqueRows[count]] + 1
                            ###############
                            # also add one to all entries of incorrectSum that have the same set of parents of i - find them
                            allOthers <- setdiff(1:(dim(mmm)[1]), uniqueRows)
                            if(length(allOthers)>1)
                            {
                                indInAllOthers <- which(colSums(!xor(mmm[uniqueRows[count],allParentsOfI],t(mmm[allOthers,allParentsOfI])))==p)
                                if(length(indInAllOthers)>0)
                                {
                                    incorrectSum[allOthers[indInAllOthers]] <- incorrectSum[allOthers[indInAllOthers]] + rep(1,length(indInAllOthers))
                                }
                            }
                            if(length(allOthers) == 1)
                            {
                                indInAllOthers <- which(sum(!xor(mmm[uniqueRows[count],allParentsOfI],t(mmm[allOthers,allParentsOfI])))==p)
                                if(length(indInAllOthers)>0)
                                {
                                    incorrectSum[allOthers[indInAllOthers]] <- incorrectSum[allOthers[indInAllOthers]] + rep(1,length(indInAllOthers))
                                }
                            }
                            #############
                            finished <- TRUE
                        }
                        
                        # if the set of parents are the same
                        if(!finished && setequal(paG,paGp))
                        {
                            finished <- TRUE
                            correctInt[i,j] <- 1
                        }
                        
                        # this part contains the difficult computations
                        if( !finished )
                        {
                            if(PathMatrix[i, j] > 0)
                            {
                                #which children are part of a causal path?
                                chiCausPath <- which( trueGraph[i,] & PathMatrix[,j])
                                #check whether in paGp there is a descendant of a "proper" child of i 
                                if(sum(PathMatrix[chiCausPath,paGp])>0)
                                {
                                    incorrectInt[i,j] <- 1
                                    incorrectSum[uniqueRows[count]] <- incorrectSum[uniqueRows[count]] + 1 
                                    ###############
                                    # also add one to all entries of incorrectSum that have the same set of parents of i - find them
                                    allOthers <- setdiff(1:(dim(mmm)[1]), uniqueRows)
                                    if(length(allOthers)>1)
                                    {
                                        indInAllOthers <- which(colSums(!xor(mmm[uniqueRows[count],allParentsOfI],t(mmm[allOthers,allParentsOfI])))==p)
                                        if(length(indInAllOthers)>0)
                                        {
                                            incorrectSum[allOthers[indInAllOthers]] <- incorrectSum[allOthers[indInAllOthers]] + rep(1,length(indInAllOthers))
                                        }
                                    }
                                    if(length(allOthers) == 1)
                                    {
                                        indInAllOthers <- which(sum(!xor(mmm[uniqueRows[count],allParentsOfI],t(mmm[allOthers,allParentsOfI])))==p)
                                        if(length(indInAllOthers)>0)
                                        {
                                            incorrectSum[allOthers[indInAllOthers]] <- incorrectSum[allOthers[indInAllOthers]] + rep(1,length(indInAllOthers))
                                        }
                                    }
                                    #############
                                    finished <- TRUE
                                }
                            }
                            
                            if(!finished)
                            {
                                #check whether all non-causal paths are blocked
                                if(reachableWOutCausalPath[j]==1)
                                {
                                    incorrectInt[i,j] <-  1
                                    incorrectSum[uniqueRows[count]] <- incorrectSum[uniqueRows[count]] + 1
                                    ###############
                                    # also add one to all entries of incorrectSum that have the same set of parents of i - find them
                                    allOthers <- setdiff(1:(dim(mmm)[1]), uniqueRows)
                                    if(length(allOthers)>1)
                                    {
                                        indInAllOthers <- which(colSums(!xor(mmm[uniqueRows[count],allParentsOfI],t(mmm[allOthers,allParentsOfI])))==p)
                                        if(length(indInAllOthers)>0)
                                        {
                                            incorrectSum[allOthers[indInAllOthers]] <- incorrectSum[allOthers[indInAllOthers]] + rep(1,length(indInAllOthers))
                                        }
                                    }
                                    if(length(allOthers) == 1)
                                    {
                                        indInAllOthers <- which(sum(!xor(mmm[uniqueRows[count],allParentsOfI],t(mmm[allOthers,allParentsOfI])))==p)
                                        if(length(indInAllOthers)>0)
                                        {
                                            incorrectSum[allOthers[indInAllOthers]] <- incorrectSum[allOthers[indInAllOthers]] + rep(1,length(indInAllOthers))
                                        }
                                    }
                                    #############
                                    
                                } else
                                {
                                    correctInt[i,j] <- 1
                                }
                            }                            
                        }
                    }
                } #for-loop over j
                count <- count + 1   
            } # while-loop count<=maxcount 
            if(!GpIsEssentialGraph)
            {
                minimumTotal <- minimumTotal + min(incorrectSum)
                maximumTotal <- maximumTotal + max(incorrectSum)
                incorrectSum <- 0
            }
            if(length(incorrectSum)>1 & output)
            {
                cat("For variable ", i, "we have more than one possible set of parents.\n")
                cat("The following vector is adding the number of incorrect interventions for the possible parent sets.\n")
                cat("(the computation is done in a clever way since some of the components may correspond to the same parent set of node ",i,").\n")
                show(incorrectSum)
                cat("Only for a new connected component this vector is set to zero again.\n\n")
            }
        } # i in conn.comp
        minimumTotal <- minimumTotal + min(incorrectSum)
        maximumTotal <- maximumTotal + max(incorrectSum)
        incorrectSum <- 0
    }
    timeTotal <- proc.time()-ptmtotal
    
    # The rest is output and return.
    if(output && p < 11)
    {
        show("These all are incorrectly predicted interventions:")
        show(incorrectInt)
        show("And these are all correctly predicted interventions:")
        show(correctInt)
    }
    ress <- list()
    ress$sid <- sum(incorrectInt)
    ress$sidUpperBound <- maximumTotal
    ress$sidLowerBound <- minimumTotal
    ress$incorrectMat <- incorrectInt
    if(output)
    {
        cat("Time needed for ... \n")
        cat("... expending the graph: ",timeExpGraph[3] ,"\n")
        cat("... computing path matrices (used for checking d-seps): ",timePathMatrix2[3] ,"\n")
        cat("... checking d-separations: ",timeAlldSep[3] ,"\n")
        cat("... in total: ",timeTotal[3] ,"\n")
        cat("number of times we ran *check all d-seps*:", numChecks,"\n")
    }
    return(ress)
}
