# Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
library(pcalg) 

GES_with_score <- function(X, X_val, lambda)
{
  result <- list()
  X <- as.matrix(X)
  X_val <- as.matrix(X_val)
  
  score <- new("GaussL0penObsScore", X, lambda=lambda)
  G <- ges(score)
  result$Adj <- as(G$essgraph, "matrix")
  n <- dim(X)[1]
  p <- dim(X)[2]
  #result$Score <- -(G$essgraph$score$global.score(G$repr) - n*p/2*log(2*pi) - n*p/2 )
  # JP 16.9.2013: I am not sure why this gives sth different than the next line 
  tmp_score <- computeScoreSEMGauss(X, as(G$repr,"matrix"))
  result$score_train <- tmp_score$score
  
  tmp_score_val <- computeScoreSEMGauss(X_val, as(G$repr,"matrix"))
  result$score_val <- tmp_score_val$score
  
  show(result$score_train)
  show(result$score_val)
  
  return(result)
}


computeScoreSEMGauss <- function(X, Adj)
{
  SigmaHat <- cov(X)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  I <- diag(rep(1,p))
  B <- diag(rep(0,p))
  eachResVar <- rep(0,p)
  
  for(node in 1:(p))
  {
    parents <- which(Adj[,node]==1)
    
    if(length(parents) == 0)
    {
      # use MLE instead of unbiased estimator of the variance
      eachResVar[node] <- var(X[,node])
    }
    else
    {
      mod <- lm(X[,node] ~ X[,parents])  
      # use MLE instead of unbiased estimator of the variance
      eachResVar[node] <- var(mod$residuals)
      B[node,parents] <- mod$coef[2:(length(parents) + 1)]
    }
  }
  
  numPars <- (sum(B!=0) + p)
  
  # sigmaHatSq <- 1/(pn) sum_k sum(State$eachResVar[k] - overallMean)^2
  # but all the means are by cosntruction of the regression zero.
  SigmaNInvHat <- diag(1/eachResVar)
  
  # the following is the negloglikelihood
  #negloglik <- n/2 * log((2*pi)^p*prod(eachResVar)) + n/2 * sum( diag( t(I-B) %*% (I-B) %*% SigmaNInvHat %*% SigmaHat) )
  negloglik <- n/2 * log((2*pi)^p*prod(eachResVar)) + n/2 * p
  penalization <- log(n)/2 * numPars
  
  if(1==0)
  {
    show("eachResVar")
    show(eachResVar)
    show("negloglik SEM")
    show(negloglik)
    show("pen SEM")
    show(penalization)
  }    
  # minimize this BIC score!
  score <- -negloglik 
  penalization <- penalization
  return(list(score = score, penalization = penalization))
}


set.seed(42);
dataset <- read.csv(file='{FOLDER}{FILE_TRAIN}', sep=",");
dataset_val <- read.csv(file='{FOLDER}{FILE_VALID}', sep=",");

result <- GES_with_score(dataset, X_val = dataset_val, lambda = {LAMBDA})
write.csv(as.matrix(result$Adj),row.names = FALSE, file = '{FOLDER}{OUTPUT}');

scores <- c(result$score_train, result$score_val)
write.csv(scores, row.names = FALSE, file = '{FOLDER}{OUTPUT2}');
