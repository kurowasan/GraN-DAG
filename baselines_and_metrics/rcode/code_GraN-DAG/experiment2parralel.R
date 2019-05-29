# Modified by the authors of "Gradient-Based Neural DAG Learning"
experiment2parralel <- function(p,n,linear,pars,expCounter, dataset.path, funcType, hparams)
{
  
  res <- rep(NA,48)
  cat("\r p =",p,", n =",n, ", exp = ", expCounter, "\r")
  
  X <- np$load(paste(dataset.path, "/data", expCounter, ".npy", sep=""))
  trueG <- np$load(paste(dataset.path, "/DAG", expCounter, ".npy", sep=""))
  trueGCPDAG <- np$load(paste(dataset.path, "/CPDAG", expCounter, ".npy", sep=""))
  
  SigmaHat <- cov(X)
  
  # Brute Force
  if(F) #p<5)
  {
    show("Brute Force")
    start_time <- Sys.time()
    resBFtmp <- BruteForce(X, "SEMIND", pars, output = FALSE)
    resBF <- resBFtmp$Adj
    end_time <- Sys.time()
  } else
  {
    cat("Brute Force is not performed since the number of variables is quite high. To avoid error messages, the output is the empty DAG.\n")
    start_time <- Sys.time()
    resBF <- matrix(0,p,p)
    end_time <- Sys.time()
  }           
  res[10] <- structIntervDist(trueG, resBF)$sid
  res[11] <- hammingDistance(trueG, resBF)
  res[12] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resBF))
  res[40] <- computeFDR(trueG, resBF)
  res[34] <- end_time - start_time
  
  
  #GDS
  if(FALSE) #(p<21)
  {
    show("GDS with tabu search")
    start_time <- Sys.time()
    resGDStmp <- GDS(X, scoreName = "SEMIND", pars, check = "checkUntilFirstMinK", startAt = "emptyGraph", kvec = c(p), output = FALSE, tabu = TRUE)
    resGDS <- resGDStmp$Adj
    end_time <- Sys.time()
  } else
  {
    start_time <- Sys.time()
    cat("GDS is not performed since the number of variables is quite high. To avoid error messages, the output is the empty DAG.\n")
    resGDS <- matrix(0,p,p)
    end_time <- Sys.time()
  }
  res[1] <- structIntervDist(trueG, resGDS)$sid
  res[2]  <- hammingDistance(trueG, resGDS)
  res[3]  <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resGDS))
  res[41] <- computeFDR(trueG, resGDS)
  res[31] <- end_time - start_time
  
  # RESIT (in the code also called ICML)
  if(p<21){
    show("RESIT")
    start_time <- Sys.time()
    resICML <- ICML(X, alpha=hparams$resit_alpha, model=pars$regr.method, indtest=indtestHsic, output = FALSE)
    end_time <- Sys.time()
  } else{
    # this method performs O(#nodes ^ 2) independence tests which take O(#samples ^ 2).
    # https://papers.nips.cc/paper/3201-a-kernel-statistical-test-of-independence.pdf
    start_time <- Sys.time()
    cat("RESIT is not performed since the number of variables is quite high. To avoid error messages, the output is the empty DAG.\n")
    resICML <- matrix(0,p,p)
    end_time <- Sys.time()
  }
  res[4] <- structIntervDist(trueG, resICML)$sid
  res[5] <- hammingDistance(trueG, resICML)
  res[6] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resICML))
  res[42] <- computeFDR(trueG, resICML)
  res[32] <- end_time - start_time
  
  # LINGAM
  show("LINGAM")
  start_time <- Sys.time()
  resLINGAM <- lingamWrap(X)$Adj
  end_time <- Sys.time()
  res[7] <- structIntervDist(trueG, resLINGAM)$sid
  res[8] <- hammingDistance(trueG, resLINGAM)
  res[9] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resLINGAM))
  res[43] <- computeFDR(trueG, resLINGAM)
  res[33] <- end_time - start_time
  
  #PC
  show("PC")
  start_time <- Sys.time()
  resPC <- pcWrap(X, 0.01, mmax = Inf)
  ssiidd <- structIntervDist(trueG, resPC)
  end_time <- Sys.time()
  res[13] <- ssiidd$sidUpperBound
  res[14] <- ssiidd$sidLowerBound
  res[15] <- hammingDistance(trueG, resPC)
  res[16] <- hammingDistance(trueGCPDAG, resPC)
  res[44] <- computeFDR(trueG, resPC)
  res[35] <- end_time - start_time
  
  #CPC
  #show("CPC")
  #start_time <- Sys.time()
  #rescPC <- cpcWrap(X, 0.01, mmax = Inf)
  #ssiidd <- structIntervDist(trueG, rescPC)
  #end_time <- Sys.time()
  #res[17] <- ssiidd$sidUpperBound
  #res[18] <- ssiidd$sidLowerBound
  #res[19] <- hammingDistance(trueG, rescPC)
  #res[20] <- hammingDistance(trueGCPDAG, rescPC)
  #res[45] <- computeFDR(trueG, rescPC)
  #res[36] <- end_time - start_time
  
  #GES
  show("GES")
  start_time <- Sys.time()
  resGES <- gesWrap(X, hparams$ges_coeff)$Adj
  ssiidd <- structIntervDist(trueG, resGES)
  end_time <- Sys.time()
  res[21] <- ssiidd$sidUpperBound
  res[22] <- ssiidd$sidLowerBound
  res[23] <- hammingDistance(trueG, resGES)
  res[24] <- hammingDistance(trueGCPDAG, resGES)
  res[46] <- computeFDR(trueG, resGES)
  res[37] <- end_time - start_time
  
  #RANDOM
  show("RANDOM")
  start_time <- Sys.time()
  resRand <- as.matrix(randomDAG(p,runif(1)))
  end_time <- Sys.time()
  res[25] <- structIntervDist(trueG, resRand)$sid
  res[26] <- hammingDistance(trueG, resRand)      
  res[27] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resRand))
  res[47] <- computeFDR(trueG, resRand)
  res[38] <- end_time - start_time
  
  #CAM
  show("CAM")
  start_time <- Sys.time()
  resCAM <- as.matrix(CAM(X, variableSel = TRUE, variableSelMethodPars = list(atLeastThatMuchSelected = hparams$pns_threshold,
atMostThatManyNeighbors = 20), pruning=TRUE, pruneMethodPars = list(cutOffPVal = hparams$pruning_threshold, numBasisFcts = 10))$Adj)
  end_time <- Sys.time()
  res[28] <- structIntervDist(trueG, resCAM)$sid
  res[29] <- hammingDistance(trueG, resCAM)
  res[30] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resCAM))
  res[48] <- computeFDR(trueG, resCAM)
  res[39] <- end_time - start_time
  
  return(res)
}
