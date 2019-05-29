# Copyright (c) 2012-2014  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

experiment2parralel <- function(p,n,pCon,linear,pars,expCounter)
{
    res <- rep(NA,27)
    cat("\r p =",p,", n =",n, ", exp = ", expCounter, "\r")

    # Generate data
    trueG <- as.matrix(randomDAG(p,pCon))
    trueGCPDAG <- dag2cpdagAdj(trueG)
    if(linear)
    {
        trueB <- randomB(trueG,0.1,2,TRUE)
        X <- sampleDataFromG(n,trueG,funcType="linear", parsFuncType=list(B=trueB,kap=1,sigmax=1,sigmay=1,output=FALSE), noiseType="normalRandomVariancesRandomExp", parsNoise=list(varMin=0.1,varMax=0.5,noiseExpVarMin=2,noiseExpVarMax=4))
        X <- as.matrix(X)
    } else
    {
        X <- sampleDataFromG(n,trueG,funcType="GAM", parsFuncType=list(kap=1,sigmax=1,sigmay=1,output=FALSE), noiseType="normalRandomVariances", parsNoise=list(noiseExp=1,varMin=1,varMax=2))
    }
    SigmaHat <- cov(X)
    
    
    # Brute Force
    if(p<5)
    {
        show("Brute Force")
        resBFtmp <- BruteForce(X, "SEMIND", pars, output = FALSE)
        resBF <- resBFtmp$Adj
    } else
    {
	cat("Brute Force is not performed since the number of variables is quite high. To avoid error messages, the output is the empty DAG.\n")
        resBF <- matrix(0,p,p)
    }           
    res[10] <- structIntervDist(trueG, resBF)$sid
    res[11] <- hammingDistance(trueG, resBF)
    res[12] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resBF)) 
    
    
    #GDS
    if(p<21)
    {
        show("GDS with tabu search")
        resGDStmp <- GDS(X, scoreName = "SEMIND", pars, check = "checkUntilFirstMinK", startAt = "emptyGraph", kvec = c(p), output = FALSE, tabu = TRUE)
        resGDS <- resGDStmp$Adj
    } else
    {
	cat("GDS is not performed since the number of variables is quite high. To avoid error messages, the output is the empty DAG.\n")
        resGDS <- matrix(0,p,p)
    }
    res[1] <- structIntervDist(trueG, resGDS)$sid
    res[2]  <- hammingDistance(trueG, resGDS)
    res[3]  <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resGDS))
   
    
    # RESIT (in the code also called ICML)
    resICML <- ICML(X, alpha=0.05, model=pars$regr.method, indtest=indtestHsic, output = FALSE)
    res[4] <- structIntervDist(trueG, resICML)$sid
    res[5] <- hammingDistance(trueG, resICML)
    res[6] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resICML)) 
    
    
    # LINGAM
    resLINGAM <- lingamWrap(X)$Adj
    res[7] <- structIntervDist(trueG, resLINGAM)$sid
    res[8] <- hammingDistance(trueG, resLINGAM)
    res[9] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resLINGAM))

    
    #PC
    resPC <- pcWrap(X, 0.01, mmax = Inf)
    ssiidd <- structIntervDist(trueG, resPC)
    res[13] <- ssiidd$sidUpperBound
    res[14] <- ssiidd$sidLowerBound
    res[15] <- hammingDistance(trueG, resPC)
    res[16] <- hammingDistance(trueGCPDAG, resPC) 
    
    
    #CPC
    rescPC <- cpcWrap(X, 0.01, mmax = Inf)
    ssiidd <- structIntervDist(trueG, rescPC)
    res[17] <- ssiidd$sidUpperBound
    res[18] <- ssiidd$sidLowerBound
    res[19] <- hammingDistance(trueG, rescPC)
    res[20] <- hammingDistance(trueGCPDAG, rescPC) 


    #GES
    resGES <- gesWrap(X)$Adj
    ssiidd <- structIntervDist(trueG, resGES)
    res[21] <- ssiidd$sidUpperBound
    res[22] <- ssiidd$sidLowerBound
    res[23] <- hammingDistance(trueG, resGES)
    res[24] <- hammingDistance(trueGCPDAG, resGES) 
    
    
    #RANDOM
    resRand <- as.matrix(randomDAG(p,runif(1)))
    res[25] <- structIntervDist(trueG, resRand)$sid
    res[26] <- hammingDistance(trueG, resRand)      
    res[27] <- hammingDistance(trueGCPDAG, dag2cpdagAdj(resRand))
    
    return(res)
}
