# Copyright (c) 2013-2013  Jonas Peters  [peters@stat.math.ethz.ch]
computeAllScoresParallel <- function(i, X, scoreName, parsScore, selMat, lookUpParents, output)
{
    
    parents <- which(selMat[,i]==TRUE)
    numParents <- length(parents)
    numParentsSets <- 2^(numParents)
    if(output)
    {
        cat("Computing score entry for ",i,"\n")
        cat("There are ", numParentsSets ,"possible parent sets for ",i,"\n")
    }
    scoresForI <- rep(NA,numParentsSets)
    if(numParents>20)
    {
        error("Parents set of size 20? This is not going to work.")
    }
    
    scoresForI[1] <- -log(var(X[,i]))
    
    if(numParentsSets > 1)
    {
        for(paset in 1:(numParentsSets-1))
        {
            if(output)
            {
                cat("\r Computing ", paset, " out of ", numParentsSets ,"possible parent sets")
            }
            parentsToCheck <- parents[which(as.integer(intToBits(paset))[1:numParents] == 1)]
            #show(parentsToCheck)
            
            if(scoreName == "SEMSEV")
            {      
                stop("This score does not work. It does not decouple.")
            } else if(scoreName == "SEMIND")
            {
                stop("NOT IMPLEMENTED YET")
            } else if(scoreName == "SEMGAM")
            {
                mod_gam <- train_gam(X[,parentsToCheck],X[,i],pars=parsScore)
                # log likelihood, without minus. We want to MAXIMIZE this score.
                scoresForI[paset+1] <- (-log(var(mod_gam$residuals)))
            } else if(scoreName == "SEMLIN")
            {
                mod_gam <- train_linear(X[,parentsToCheck],X[,i])
                # log likelihood, without minus. We want to MAXIMIZE this score.
                scoresForI[paset+1] <- (-log(var(mod_gam$residuals)))
            } else if(scoreName == "SEMGP")
            {
                mod_gp <- train_gp(X[,parentsToCheck],X[,i])
                scoresForI[paset+1] <- (-log(var(mod_gp$residuals)))
            } else
            {
                stop("I do not know this score function.")
            }
        }
    }
    
    if(output)
    {
        cat("\n Done with ", i ,"\n")
    }
    
    return(scoresForI)
}


computeAllScores <-function(X, scoreName, output, numCores, selMat, parsScore, lookUpParents)
{
    p <- dim(X)[2]
    n <- dim(X)[1]
    
    allNodes <- as.list(1:p)
    
    computeAllScoresParallelNew <- function(i){return(computeAllScoresParallel(i, X = X, scoreName = scoreName, parsScore = parsScore, selMat = selMat, lookUpParents = lookUpParents, output = output))}
    
    if(numCores == 1)
    {
        allScores <- lapply(allNodes,computeAllScoresParallelNew)
    } else
    {
        allScores <- mclapply(allNodes,computeAllScoresParallelNew, mc.cores = numCores)
    }
    
    return(allScores)
}