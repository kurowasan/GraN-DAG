kPaGreedy <- function(X, scoreName, k, output = FALSE) 
{
    timeCycle <- 0
    timeScoreMat <- 0
    ptm <- proc.time()
    temp <- computeScoreMat(X, scoreName, k, output) 
    timeScoreMat <- timeScoreMat + proc.time() - ptm
    scoreMat <- temp$scoreMat
    rowParents <- temp$rowParents
    p <- dim(scoreMat)[2]
    Adj <- matrix(0,p,p)
    counter <- 0
    score <- 0
    parentsNotFoundYet <- 1:p
    while( counter < (p-k) ) 
    {
        
        ix_max_scoreMat <- arrayInd(which.max(scoreMat), dim(scoreMat))
        
        #show(Adj)
        #show(scoreMat)
        #cat(rowParents[ix_max_scoreMat[1],]," on ", ix_max_scoreMat[2],"\n")
        #readline()
        
        Pa <- rowParents[ix_max_scoreMat[1],]
        Chi <- rep(ix_max_scoreMat[2], length(Pa))
        ix_max_Adj <- t(rbind(Pa, Chi))
        Adj[ix_max_Adj] <- 1
        
        ptm <- proc.time()[3]
        cCA <- containsCycle(Adj)
        timeCycle <- timeCycle + proc.time() - ptm
        
        if (cCA == TRUE)
        {
            Adj[ix_max_Adj] <- 0
            scoreMat[ix_max_scoreMat[1],ix_max_scoreMat[2]] <- -Inf
        } else {
            score <- score + scoreMat[ix_max_scoreMat[1],ix_max_scoreMat[2]]
            scoreMat[, ix_max_scoreMat[2]] <- -Inf
            counter <- counter + 1
            parentsNotFoundYet <- setdiff(parentsNotFoundYet, ix_max_scoreMat[2])
        }
    }
    
    if(k > 1) 
    {
        if(output)
        {
            "finding DAG for remaining nodes..."
        }
        ttmp <- BruteForceFullDags(X[,parentsNotFoundYet], scoreName, pars = list(), output)
        score <- score - ttmp$score
        remainingDAG <- ttmp$Adj 
        if(sum(Adj[parentsNotFoundYet,parentsNotFoundYet])> 0)
        {
            stop("BUG ATTACK!")
        }
        Adj[parentsNotFoundYet,parentsNotFoundYet] <- remainingDAG
    }
    show(timeCycle)
    show(timeScoreMat)
    result <- list(Adj = Adj, Score = score)
    return(result)
}
