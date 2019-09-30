# MIT License
#
# Copyright (c) 2018 Diviyan Kalainathan
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

library(CAM)

CAM_with_score <-
  function(X, scoreName = "SEMGAM", 
           parsScore = list(numBasisFcts=10), 
           numCores = 1, 
           maxNumParents = min(dim(X)[2] - 1, round(dim(X)[1]/20)),
           output = FALSE, 
           variableSel = FALSE, 
           variableSelMethod = selGamBoost, 
           variableSelMethodPars = list(atLeastThatMuchSelected = 0.02, atMostThatManyNeighbors = 10),
           is_pruning = FALSE, 
           pruneMethod = selGam, 
           pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts=10),
           intervData = FALSE, 
           intervMat = NA,
           X_val=NULL) 
  {
    require(CAM)
    
    if(output)
    {
      cat("number of cores:", numCores, "\n")
    }
    # We record the time consumption. They are shown if output == TRUE
    timeCycle <- 0
    timeUpdate <- 0
    timeScoreMat <- 0
    timeSel <- 0
    timePrune <- 0
    timeMax <- 0
    
    # we record how the score develops 
    scoreVec <- integer(0)
    scoreVec_val <- integer(0)
    # and which edges are added
    edgeList <- integer(0)
    
    # this counter is only used if output = TRUE
    counterUpdate <- 0
    p <- dim(X)[2]
    
    
    ####
    # STEP 1: variable selection
    ####
    # A matrix selMat is constructed. Entry (i,j) being one means that i is a possible parent of j.
    if(variableSel)
    {
      ptm <- proc.time()[3]
      if(intervData)
      {
        X2 <- X[rowSums(intervMat) == 0,]
        if(output)
          cat("The preliminary neighbourhood selection is done with the observational data only.\n")
      } else
      {
        X2 <- X
      }
      if(numCores == 1)
      {
        selMat <- mapply(variableSelMethod,MoreArgs = list(X = X2, pars = variableSelMethodPars, output = output),1:p)
      } else
      {
        selMat <- mcmapply(variableSelMethod,MoreArgs = list(X = X2, pars = variableSelMethodPars, output = output),1:p, mc.cores = numCores)
      }
      # The next line includes j as a possible parent of i if i is considered a possible parent of j
      # selMat <- selMat | t(selMat)
      cou <- 0
      for(jk in 1:p)
      {
        cou <- cou + 2^{sum(selMat[,jk])}
      }
      if(output)
      {
        cat("Instead of p2^(p-1) -Sillander- ",p*2^(p-1) ," we have ", cou, "\n")
        cat("Greedy, on the other hand, is computing ",sum(selMat) ," entries. \n")
      }
      timeSel <- timeSel + proc.time()[3] - ptm
    } else
    {
      selMat <- matrix(TRUE, p,p)
    }
    if(variableSel & output)
    {
      if(output)
      {
        if(p<30)
        {
          cat("This is the matrix of possible parents after the first step.\n")
          show(selMat)
        }
        cat("Object size of selmat: ", object.size(selMat), "\n")
      }
    }
    
    
    ####
    # STEP 2: Include Edges
    ####
    # compute score matrix 
    ptm <- proc.time()[3]
    computeScoreMatTmp <- computeScoreMat_with_score(X, X_val, scoreName=scoreName, numParents = 1, numCores = numCores, output = output, selMat = selMat, parsScore = parsScore, intervMat = intervMat, intervData = intervData)
    original_score_mat <- computeScoreMatTmp$scoreMatOriginal
    original_score_mat_val <- computeScoreMatTmp$scoreMat_val
    
    
    timeScoreMat <- timeScoreMat + proc.time()[3] - ptm
    if(output)
    {
      cat("Object size of computeScoreMatTmp: ", object.size(computeScoreMatTmp), "\n" )
    }
    # We need the pathMatrix (entry (i,j) being one means that there is a directed path from i to j) in order to keep track of possible cycles.
    pathMatrix <- matrix(0,p,p)
    diag(pathMatrix) <- rep(1,p)
    Adj <- as(matrix(0,p,p), "sparseMatrix")
    scoreNodes <- computeScoreMatTmp$scoreEmptyNodes
    scoreNodes_val <- computeScoreMatTmp$scoreEmptyNodes_val
    i <- 0
    
    # Greedily adding edges
    while(sum(computeScoreMatTmp$scoreMat!=-Inf) > 0)
    {
      # Find the best edge
      ptm <- proc.time()[3]
      ix_max <- arrayInd(which.max(computeScoreMatTmp$scoreMat), dim(computeScoreMatTmp$scoreMat))
      ix_max_backward <- matrix(c(ix_max[2],ix_max[1]),1,2)
      if(i == 0){
        order_ix_max <- data.frame(t(ix_max))
      }else{
        order_ix_max <- cbind(order_ix_max, t(ix_max))
      }
      
      timeMax <- timeMax + proc.time()[3] - ptm
      Adj[ix_max] <- 1
      scoreNodes[ix_max[2]] <- scoreNodes[ix_max[2]] + computeScoreMatTmp$scoreMat[ix_max]
      scoreNodes_val[ix_max[2]] <- scoreNodes_val[ix_max[2]] + computeScoreMatTmp$scoreMat_val[ix_max]

      if(output)
      {
        cat("\n Included edge (from, to) ", ix_max, "\n")
      }
      
      # Do not include the same edge twice.
      computeScoreMatTmp$scoreMat[ix_max] <- -Inf
      computeScoreMatTmp$scoreMat_val[ix_max] <- -Inf
      
      
      # Avoid cycles
      ptm <- proc.time()[3]
      pathMatrix[ix_max[1],ix_max[2]] <- 1
      DescOfNewChild <- which(pathMatrix[ix_max[2],]==1)
      AncOfNewParent <- which(pathMatrix[,ix_max[1]]==1)
      pathMatrix[AncOfNewParent,DescOfNewChild] <- 1
      computeScoreMatTmp$scoreMat[t(pathMatrix) == 1] <- -Inf 
      computeScoreMatTmp$scoreMat[ix_max[2],ix_max[1]] <- -Inf
      computeScoreMatTmp$scoreMat_val[t(pathMatrix) == 1] <- -Inf 
      computeScoreMatTmp$scoreMat_val[ix_max[2],ix_max[1]] <- -Inf
      
      timeCycle <- timeCycle + proc.time()[3] - ptm
      
      # Record the score of the current graph
      scoreVec <- c(scoreVec, sum(scoreNodes))
      scoreVec_val <- c(scoreVec_val, sum(scoreNodes_val))
      # Record which edge has been added
      edgeList <- rbind(edgeList, ix_max, deparse.level=0)
      
      # Update the score of column j
      ptm <- proc.time()[3]
      updates <- updateScoreMat_with_score(computeScoreMatTmp$scoreMat, original_score_mat,  computeScoreMatTmp$scoreMat_val, X, X_val, scoreName = scoreName, ix_max[1], ix_max[2], scoreNodes, scoreNodes_val, Adj, numCores=numCores, output = output, maxNumParents = maxNumParents, parsScore = parsScore, intervMat = intervMat, intervData = intervData)
      computeScoreMatTmp$scoreMat <- updates$scoreMat
      original_score_mat <- updates$scoreMat_original
      computeScoreMatTmp$scoreMat_val <- updates$scoreMat_val

      timeUpdate <- timeUpdate + proc.time()[3] - ptm
      
      counterUpdate <- counterUpdate + 1
      i <- i + 1
      
    }
    
    # print("Score before pruning:")
    # print(tail(scoreVec, n=1))
    # print(tail(scoreVec_val, n=1))

    Adj_prepruned <- as.matrix(Adj)
    

    ####
    # STEP 3: Prune the DAG
    ####
    if(is_pruning)
    {
      if(intervData)
      {
        X2 <- X[rowSums(intervMat) == 0,]
        cat("The preliminary neighbourhood selection is done with the observational data only.\n")
      } else
      {
        X2 <- X
      }
      if(output)
      {
        cat("\n Performing pruning ... \n ")
      }
      ptm <- proc.time()[3]
      Adj <- CAM:::pruning(X=X2,G=Adj,pruneMethod = pruneMethod, pruneMethodPars = pruneMethodPars, output=output)  
      
      Adj_tmp <- as(matrix(0,p,p), "sparseMatrix")
      computeScoreMatTmp <- computeScoreMat_with_score(X, X_val, scoreName=scoreName, numParents = 1, numCores = numCores, output = output, selMat = selMat, parsScore = parsScore, intervMat = intervMat, intervData = intervData)
      original_score_mat <- computeScoreMatTmp$scoreMatOriginal
      original_score_mat_val <- computeScoreMatTmp$scoreMat_val
      timePrune <- timePrune + proc.time()[3] - ptm    
      

      #Recompute the score after pruning
      scoreNodes <- computeScoreMatTmp$scoreEmptyNodes
      scoreNodes_val <- computeScoreMatTmp$scoreEmptyNodes_val
      scoreVec <- c()
      scoreVec_val <- c()
      
      for(i in 1:dim(order_ix_max)[2]){
        ix_max <- order_ix_max[,i]
        
        if(Adj[ix_max[1], ix_max[2]] == 1){
            Adj_tmp[ix_max[1], ix_max[2]] <- 1
            ix_max_backward <- matrix(c(ix_max[2],ix_max[1]),1,2)

            scoreNodes[ix_max[2]] <- scoreNodes[ix_max[2]] + computeScoreMatTmp$scoreMat[ix_max[1], ix_max[2]]
            scoreNodes_val[ix_max[2]] <- scoreNodes_val[ix_max[2]] + computeScoreMatTmp$scoreMat_val[ix_max[1], ix_max[2]]
            
            # Record the score of the current graph
            scoreVec <- c(scoreVec, sum(scoreNodes))
            scoreVec_val <- c(scoreVec_val, sum(scoreNodes_val))
            
            # Update the score of column j
            updates <- updateScoreMat_with_score(computeScoreMatTmp$scoreMat, original_score_mat,  computeScoreMatTmp$scoreMat_val, X, X_val, scoreName = scoreName, ix_max[1], ix_max[2], scoreNodes, scoreNodes_val, Adj_tmp, numCores=numCores, output = output, maxNumParents = maxNumParents, parsScore = parsScore, intervMat = intervMat, intervData = intervData)
            computeScoreMatTmp$scoreMat <- updates$scoreMat
            original_score_mat <- updates$scoreMat_original
            computeScoreMatTmp$scoreMat_val <- updates$scoreMat_val
        }
      }
    }
    
    # print("Score after pruning:")
    # print(tail(scoreVec, n=1))
    # print(tail(scoreVec_val, n=1))
    
    final_score <- tail(scoreVec, n=1)
    final_score_val <- tail(scoreVec_val, n=1)
    
    show(paste("number of edges: ", sum(Adj), sep=""))
    
    
    ####
    # Output and return
    ####
    timeTotal <- timeSel + timeScoreMat + timeCycle + timeUpdate + timeMax + timePrune
    if(output)
    {
      cat("amount of time for variable selection:",timeSel,"\n")
      cat("amount of time computing the initial scoreMat:",timeScoreMat,"\n")
      cat("amount of time checking for cycles:",timeCycle,"\n")
      cat("amount of time computing updates for the scoreMat:",timeUpdate,", doing",counterUpdate,"updates.\n")
      cat("amount of time for pruning:",timePrune,"\n")
      cat("amount of time for finding maximum:",timeMax,"\n")
      cat("amount of time in total:",timeTotal,"\n")
    }
    
    result <- list(Adj = Adj, Score = sum(scoreNodes), timesVec = c(timeSel, timeScoreMat, timeCycle, timeUpdate, timePrune, timeMax, timeTotal), scoreVec = scoreVec, edgeList = edgeList, final_score = final_score, final_score_val = final_score_val)
    return(result)  
  }


computeScoreMat_with_score <-
  function(X, X_val, scoreName, numParents, output, numCores, selMat, parsScore, intervMat, intervData)
  {
    
    # numParents indicates how many parents we consider. If numParents = 1 (default), then the 
    # score matrix is of dimension (p-1) x p. If numParents = 2, then the  
    # score matrix is of dimension (p-1)(p-2) x p and so on...
    #
    # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
    # it should therefore be positive.
    
    p <- dim(X)[2]
    n <- dim(X)[1]
    rowParents <- t(combn(p,numParents))
    
    tt <- expand.grid(1:dim(rowParents)[1], 1:p)
    allNode2 <- tt[,2]
    allI <- tt[,1]
    if(numCores == 1)
    {
      scoreMat <- mapply(computeScoreMatParallel_with_score, MoreArgs = list(rowParents = rowParents, selMat = selMat, scoreName = scoreName, X = X, X_val = X_val, output = output, parsScore = parsScore, intervMat = intervMat, intervData = intervData), node2 = allNode2, i = allI)
    } else
    {
      scoreMat <- mcmapply(computeScoreMatParallel_with_score, MoreArgs = list(rowParents = rowParents, selMat = selMat, scoreName = scoreName, X = X, X_val = X_val, output = output, parsScore = parsScore, intervMat = intervMat, intervData = intervData), node2 = allNode2, i = allI, mc.cores = numCores)
    }
    
    
    print(dim(scoreMat))
    scoreMat_train <- matrix(unlist(scoreMat[1,]),dim(rowParents)[1],p)
    scoreMat_train_original <-matrix(unlist(scoreMat[1,]),dim(rowParents)[1],p)
    scoreMat_val <- matrix(unlist(scoreMat[2,]),dim(rowParents)[1],p)
    
    # initScore[i] equals the variance of variable j. 
    initScore <- rep(NA,p)
    for(i in 1:p)
    {
      if(intervData)
      {
        X2 <- X[!intervMat[,i],]
      } else
      {
        X2 <- X
      }
      vartmp <- var(X2[,i])
      initScore[i] <- -log(vartmp)
      # scoreMat[i,j] equals the GAIN in score if we consider i being a parent of j. 
      scoreMat_train[,i] <- scoreMat_train[,i] - initScore[i]
    }

    initScore_val <- rep(NA,p)
    for(i in 1:p)
    {
      X2 <- X_val
      vartmp_val <- var(X2[,i])
      initScore_val[i] <- -log(vartmp_val)
      scoreMat_val[,i] <- scoreMat_val[,i] - initScore_val[i]
    }

    return(list(scoreMat = scoreMat_train, scoreMatOriginal = scoreMat_train_original, scoreMat_val = scoreMat_val, rowParents = rowParents, scoreEmptyNodes = initScore, scoreEmptyNodes_val = initScore_val))
  }


computeScoreMatParallel_with_score <-
  function(rowParents, scoreName, X, X_val, selMat, output, node2, i, parsScore, intervMat, intervData, update=FALSE)
  {
    #the i-th row of rowParents contains possible parents of node2 (we call them "parentsToCheck") 
    parentsToCheck <- rowParents[i,]
    if(output)
    {
      cat("\r compute score entry for regressing",node2,"on",parentsToCheck,"                  \r")
    }
    if(intervData)
    {
      X2 <- X[!intervMat[,node2],]
    } else
    {
      X2 <- X
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
        mod_gam <- train_gam_with_score(X2[,parentsToCheck],X2[,node2],pars=parsScore)
        score <- (-log(var(mod_gam$residuals)))

        nb_var <- count(mod_gam$formula, "+") + 1
        data_frame <- data.frame(X_val[,parentsToCheck])
        var_name <- c()
        for (v in 1:nb_var){
          var_name <- c(var_name, paste("var", v+1, sep="")) 
        }
        colnames(data_frame) <- c(var_name)
        
        val <- predict.gam(mod_gam$model, data_frame)
        residual <- X_val[,node2] - val
        score_val <- (-log(var(residual)))
      } else if(scoreName == "SEMLIN")
      {
        mod_gam <- train_linear(X2[,parentsToCheck],X2[,node2])
        score <- (-log(var(mod_gam$residuals)))
      } else if(scoreName == "SEMGP")
      {
        mod_gp <- train_gp(X2[,parentsToCheck],X2[,node2])
        score <- (-log(var(mod_gp$residuals)))
      } else
      {
        stop("I do not know this score function.")
      }
    } else
    {
      score <- (-Inf)
      score_val <- (-Inf)
    }

      
    return(list(score_train = score, score_val = score_val))
  }

count = function(haystack, needle)
{v = attr(gregexpr(needle, haystack, fixed = T)[[1]], "match.length")
if (identical(v, -1L)) 0 else length(v)}

train_gam_with_score <-
  function(X,y,pars = list(numBasisFcts = 10))
  {

    if(!("numBasisFcts" %in% names(pars) ))
    { 
      pars$numBasisFcts = 10
    }
    p <- dim(as.matrix(X))
    if(p[1]/p[2] < 3*pars$numBasisFcts)
    {
      pars$numBasisFcts <- ceiling(p[1]/(3*p[2]))
      cat("changed number of basis functions to    ", pars$numBasisFcts, "    in order to have enough samples per basis function\n")
    }
    dat <- data.frame(as.matrix(y),as.matrix(X))
    coln <- rep("null",p[2]+1)
    for(i in 1:(p[2]+1))
    {
      coln[i] <- paste("var",i,sep="")
    }
    colnames(dat) <- coln
    labs<-"var1 ~ "
    if(p[2] > 1)
    {
      for(i in 2:p[2])
      {
        labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,") + ",sep="")
      }
    }
    labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,")",sep="")
    mod_gam <- FALSE
    try(mod_gam <- gam(formula=formula(labs), data=dat),silent = TRUE)
    if(typeof(mod_gam) == "logical")
    {
      cat("There was some error with gam. The smoothing parameter is set to zero.\n")
      labs<-"var1 ~ "
      if(p[2] > 1)
      {
        for(i in 2:p[2])
        {
          labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,",sp=0) + ",sep="")
        }
      }
      labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,",sp=0)",sep="")
      mod_gam <- gam(formula=formula(labs), data=dat)
    }
    result <- list()
    result$Yfit <- as.matrix(mod_gam$fitted.values)
    result$residuals <- as.matrix(mod_gam$residuals)
    result$model <- mod_gam 
    result$df <- mod_gam$df.residual     
    result$edf <- mod_gam$edf     
    result$edf1 <- mod_gam$edf1
    result$formula <- labs
    
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
    return(result)
  }

updateScoreMat_with_score <-
  function(scoreMat, scoreMat_train_original, scoreMat_val, X, X_val, scoreName, i, j, scoreNodes, scoreNodes_val, Adj, output, numCores, maxNumParents, parsScore, intervMat, intervData)
    # new edge: from i to j
  {
    p <- dim(X)[2]
    existingParOfJ <- which(Adj[,j] == 1)
    notAllowedParOfJ <- setdiff(which(scoreMat[,j] == -Inf), c(existingParOfJ,j))
    
    # if there is something left that we need to update
    if(length(existingParOfJ) + length(notAllowedParOfJ) < p-1)
    {
      # update column for j
      rowParents <- matrix(c(existingParOfJ,NA), p, length(existingParOfJ)+1, byrow = TRUE)
      rowParents[,length(existingParOfJ)+1] <- 1:p
      toUpdate <- setdiff(1:p,c(j,existingParOfJ,notAllowedParOfJ))

      
      if(length(existingParOfJ)< maxNumParents)
      {
        if(numCores == 1)
        {
          scoreUpdate <- mapply(computeScoreMatParallel_with_score,MoreArgs = list(rowParents = rowParents, selMat = matrix(TRUE,p,p), scoreName = scoreName, X = X, X_val = X_val, output = output, node2 = j, parsScore = parsScore, intervMat = intervMat, intervData = intervData, update=TRUE), i = toUpdate)
        } else
        {
          scoreUpdate <- mcmapply(computeScoreMatParallel_with_score,MoreArgs = list(rowParents = rowParents, selMat = matrix(TRUE,p,p), scoreName = scoreName, X = X, X_val = X_val, output = output, node2 = j, parsScore = parsScore, intervMat = intervMat, intervData = intervData, update=TRUE), i = toUpdate, mc.cores = numCores)
        }
        score_train <- scoreUpdate['score_train',]
        score_val <- scoreUpdate['score_val',]
        
        scoreUpdate$score <- unlist(score_train)
        scoreUpdate$score_val <- unlist(score_val)
      } else
      {
        scoreUpdate$score <- -Inf - scoreNodes[j]
        scoreUpdate$score_val <- -Inf - scoreNodes_val[j]
      }
      
      scoreMat_train_original[toUpdate,j] <- scoreUpdate$score
      scoreMat[toUpdate,j] <- scoreUpdate$score - scoreNodes[j]
      
      scoreMat_val[toUpdate,j] <- scoreUpdate$score_val - scoreNodes_val[j]
    }    
    return(list(scoreMat = scoreMat, scoreMat_original = scoreMat_train_original, scoreMat_val = scoreMat_val))
  }

set.seed(42)
dataset <- read.csv(file='{FOLDER}{FILE_TRAIN}', sep=",");
dataset_val <- read.csv(file='{FOLDER}{FILE_VALID}', sep=",");

estDAG <- CAM_with_score(dataset, X_val = dataset_val, scoreName = "{SCORE}", numCores = {NJOBS}, output = {VERBOSE},
		      variableSel = {VARSEL}, variableSelMethod = {SELMETHOD}, is_pruning = {PRUNING},
		      pruneMethod = {PRUNMETHOD}, pruneMethodPars = list(cutOffPVal = {CUTOFF}))
write.csv(as.matrix(estDAG$Adj),row.names = FALSE, file = '{FOLDER}{OUTPUT}');

scores <- c(estDAG$final_score, estDAG$final_score_val)
write.csv(scores, row.names = FALSE, file = '{FOLDER}{OUTPUT2}');
