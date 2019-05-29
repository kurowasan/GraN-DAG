computeScoreMatParallel <- function(rowParents, scoreName, X, selMat, output, node2, i, parsScore)
{
    #the i-th row of rowParents contains possible parents of node2 (we call them "parentsToCheck") 
    parentsToCheck <- rowParents[i,]
    if(output)
    {
        cat("\r compute score entry for regressing",node2,"on",parentsToCheck,"                  \r")
    }
    if(!(node2 %in% parentsToCheck) && (prod(selMat[parentsToCheck,node2]) == TRUE))
    {
        if(scoreName == "SEMSEV")
        {      
            stop("This score does not work. It does not decouple.")
        } else if(scoreName == "SEMIND")
        {
            stop("NOT IMPLEMENTED YET")
        } else if(scoreName == "SEMGAM")
        {
            mod_gam <- train_gam(X[,parentsToCheck],X[,node2],pars=parsScore)
            # log likelihood, without minus. We want to MAXIMIZE this score.
            #    scoreMat[i,node2] <- -log(var(mod_gam$residuals))
            score <- (-log(var(mod_gam$residuals)))
        } else if(scoreName == "SEMLIN")
        {
            mod_gam <- train_linear(X[,parentsToCheck],X[,node2])
            # log likelihood, without minus. We want to MAXIMIZE this score.
            #    scoreMat[i,node2] <- -log(var(mod_gam$residuals))
            score <- (-log(var(mod_gam$residuals)))
        } else if(scoreName == "SEMGP")
        {
            mod_gp <- train_gp(X[,parentsToCheck],X[,node2])
            #    scoreMat[i,node2] <- -log(var(mod_gam$residuals))
            score <- (-log(var(mod_gp$residuals)))
        } else
        {
            stop("I do not know this score function.")
        }
    } else
    {
        score <- (-Inf)
    }
    return(score)
}


computeScoreMat <- function(X, scoreName, numParents, output, numCores, selMat, parsScore)
{
    # numParents indicates how many parents we consider. If numParents = 1 (default), then the 
    # score matrix is of dimension (p-1) x p. If numParents = 2, then the  
    # score matrix is of dimension (p-1)(p-2) x p and so on...
    #
    # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
    # it should therefore be positive.
    #
    p <- dim(X)[2]
    n <- dim(X)[1]
    rowParents <- t(combn(p,numParents))
    
    tt <- expand.grid(1:dim(rowParents)[1], 1:p)
    allNode2 <- tt[,2]
    allI <- tt[,1]
    if(numCores == 1)
    {
        scoreMat <- mapply(computeScoreMatParallel,MoreArgs = list(rowParents = rowParents, selMat = selMat, scoreName = scoreName, X = X, output = output, parsScore = parsScore),node2 = allNode2, i = allI)
    } else
    {
        scoreMat <- mcmapply(computeScoreMatParallel,MoreArgs = list(rowParents = rowParents, selMat = selMat, scoreName = scoreName, X = X, output = output, parsScore = parsScore),node2 = allNode2, i = allI, mc.cores = numCores)
    }
    
    scoreMat <- matrix(scoreMat,dim(rowParents)[1],p)
    # initScore[i] equals the variance of variable j. 
    initScore <- rep(NA,p)
    for(i in 1:p)
    {
        vartmp <- var(X[,i])
        initScore[i] <- -log(vartmp)
        # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
        scoreMat[,i] <- scoreMat[,i] - initScore[i]
    }
    return(list(scoreMat = scoreMat, rowParents = rowParents, scoreEmptyNodes = initScore))
}



updateScoreMat <- function(scoreMat, X, scoreName, i, j, scoreNodes, Adj, output, numCores, maxNumParents, parsScore)
    # new edge: from i to j
{
    p <- dim(X)[2]
    existingParOfJ <- which(Adj[,j] == 1)
    notAllowedParOfJ <- setdiff(which(scoreMat[,j] == -Inf),c(existingParOfJ,j))
    
    # if there is something left that we need to update
    if(length(existingParOfJ) + length(notAllowedParOfJ) < p-1)
    {
        # update column for j
        rowParents <- matrix(c(existingParOfJ,NA),p,length(existingParOfJ)+1, byrow = TRUE)
        rowParents[,length(existingParOfJ)+1] <- 1:p
        toUpdate <- setdiff(1:p,c(j,existingParOfJ,notAllowedParOfJ))
        if(length(existingParOfJ)< maxNumParents)
        {
            if(numCores == 1)
            {
                scoreUpdate <- mapply(computeScoreMatParallel,MoreArgs = list(rowParents = rowParents, selMat = matrix(TRUE,p,p), scoreName = scoreName, X = X, output = output, node2 = j, parsScore = parsScore), i = toUpdate)
            } else
            {
                scoreUpdate <- mcmapply(computeScoreMatParallel,MoreArgs = list(rowParents = rowParents, selMat = matrix(TRUE,p,p), scoreName = scoreName, X = X, output = output, node2 = j, parsScore = parsScore), i = toUpdate, mc.cores = numCores)
            }
        } else
        {
            scoreUpdate <- -Inf
        }
        if(output)
        {
            if(sum(scoreMat[toUpdate,j] > (-log(0.8)+(scoreUpdate - scoreNodes[j]))) > 0)
            {
                cat('\n######################\n')
                cat('VIOLATION OF SUBMODULARITY: ',toUpdate[which(scoreMat[toUpdate,j] > (-log(0.8) + (scoreUpdate - scoreNodes[j])))] ,'to node',j,' \n')
                cat('######################\n') 
                show(toUpdate)
                show(scoreMat[toUpdate,j])
                show(scoreUpdate - scoreNodes[j])
            }
        }
        scoreMat[toUpdate,j] <- scoreUpdate - scoreNodes[j]
    }    
    return(scoreMat)
}

