# Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
library(mgcv)
library(MASS)
library(parallel)
source("../../util/computeGaussKernel.R", chdir = TRUE)
source("../../util_DAGs/computeCausOrder.R", chdir = TRUE)
source("./experiment2parralel.R", chdir = TRUE)
source("../../startups/startupGES.R", chdir = TRUE)
source("../../util_DAGs/randomB.R")
source("../../util_DAGs/randomDAG.R")
source("../../util_DAGs/dag2cpdagAdj.R")
source("../../util_DAGs/sampleDataFromG.R")
source("../../startups/startupSHD.R", chdir = TRUE)
source("../../startups/startupSID.R", chdir = TRUE)
source("../../startups/startupLINGAM.R", chdir = TRUE)
source("../../startups/startupICML.R", chdir = TRUE)
source("../../startups/startupBF.R", chdir = TRUE)
source("../../startups/startupGDS.R", chdir = TRUE)
source("../../startups/startupPC.R", chdir = TRUE)
source("../../startups/startupScoreSEMIND.R", chdir = TRUE)
pars <- list(regr.method = train_gam, regr.pars = list(), indtest.method = indtestHsic, indtest.pars = list())

nVec <- c(100,500)
pVec <- c(4)
numExp <- 100
numCores <-1


timeBF <- 0
timeGDS <- 0
timeICML <- 0
timeLINGAM <- 0
timeGES <- 0
timePC <- 0
timeCPC <- 0
timeRand <- 0


hdlingamDAG <- matrix(-1, length(nVec), length(pVec))
hdlingamDAGsd <- matrix(-1, length(nVec), length(pVec))
hdicmlDAG <- matrix(-1, length(nVec), length(pVec))
hdicmlDAGsd <- matrix(-1, length(nVec), length(pVec))
hdgdsDAG <- matrix(-1, length(nVec), length(pVec))
hdgdsDAGsd <- matrix(-1, length(nVec), length(pVec))
hdbfDAG <- matrix(-1, length(nVec), length(pVec))
hdbfDAGsd <- matrix(-1, length(nVec), length(pVec))
hdcpcDAG <- matrix(-1, length(nVec), length(pVec))
hdcpcDAGsd <- matrix(-1, length(nVec), length(pVec))
hdpcDAG <- matrix(-1, length(nVec), length(pVec))
hdpcDAGsd <- matrix(-1, length(nVec), length(pVec))
hdgesDAG <- matrix(-1, length(nVec), length(pVec))
hdgesDAGsd <- matrix(-1, length(nVec), length(pVec))
hdrandDAG <- matrix(-1, length(nVec), length(pVec))
hdrandDAGsd <- matrix(-1, length(nVec), length(pVec))

hdlingamCPDAG <- matrix(-1, length(nVec), length(pVec))
hdlingamCPDAGsd <- matrix(-1, length(nVec), length(pVec))
hdicmlCPDAG <- matrix(-1, length(nVec), length(pVec))
hdicmlCPDAGsd <- matrix(-1, length(nVec), length(pVec))
hdgdsCPDAG <- matrix(-1, length(nVec), length(pVec))
hdgdsCPDAGsd <- matrix(-1, length(nVec), length(pVec))
hdbfCPDAG <- matrix(-1, length(nVec), length(pVec))
hdbfCPDAGsd <- matrix(-1, length(nVec), length(pVec))
hdcpcCPDAG <- matrix(-1, length(nVec), length(pVec))
hdcpcCPDAGsd <- matrix(-1, length(nVec), length(pVec))
hdpcCPDAG <- matrix(-1, length(nVec), length(pVec))
hdpcCPDAGsd <- matrix(-1, length(nVec), length(pVec))
hdgesCPDAG <- matrix(-1, length(nVec), length(pVec))
hdgesCPDAGsd <- matrix(-1, length(nVec), length(pVec))
hdrandCPDAG <- matrix(-1, length(nVec), length(pVec))
hdrandCPDAGsd <- matrix(-1, length(nVec), length(pVec))

sidlingamDAG <- matrix(-1, length(nVec), length(pVec))
sidlingamDAGsd <- matrix(-1, length(nVec), length(pVec))
sidicmlDAG <- matrix(-1, length(nVec), length(pVec))
sidicmlDAGsd <- matrix(-1, length(nVec), length(pVec))
sidgdsDAG <- matrix(-1, length(nVec), length(pVec))
sidgdsDAGsd <- matrix(-1, length(nVec), length(pVec))
sidbfDAG <- matrix(-1, length(nVec), length(pVec))
sidbfDAGsd <- matrix(-1, length(nVec), length(pVec))
sidcpcDAGU <- matrix(-1, length(nVec), length(pVec))
sidcpcDAGUsd <- matrix(-1, length(nVec), length(pVec))
sidcpcDAGL <- matrix(-1, length(nVec), length(pVec))
sidcpcDAGLsd <- matrix(-1, length(nVec), length(pVec))
sidpcDAGU <- matrix(-1, length(nVec), length(pVec))
sidpcDAGUsd <- matrix(-1, length(nVec), length(pVec))
sidpcDAGL <- matrix(-1, length(nVec), length(pVec))
sidpcDAGLsd <- matrix(-1, length(nVec), length(pVec))
sidgesDAGU <- matrix(-1, length(nVec), length(pVec))
sidgesDAGUsd <- matrix(-1, length(nVec), length(pVec))
sidgesDAGL <- matrix(-1, length(nVec), length(pVec))
sidgesDAGLsd <- matrix(-1, length(nVec), length(pVec))
sidrandDAG <- matrix(-1, length(nVec), length(pVec))
sidrandDAGsd <- matrix(-1, length(nVec), length(pVec))


for(nCounter in 1:length(nVec))
{
    n <- nVec[nCounter]
    for(pCounter in 1:length(pVec))
    {
        p <- pVec[pCounter]
        #pCon <- 0.3 #non-sparse
        #pCon <- 1.5/(p-1) #sparse => 0.75p edges
        pCon <- 2/(p-1) #sparse => p edges
        
        noiseVar <- rep(1,p)
        
        if(numCores > 1)
        {
            respar <- mcmapply(experiment2parralel,MoreArgs=list(p=p,n=n,pCon=pCon,linear=FALSE, pars = pars),1:numExp, mc.cores = numCores)
        } else
        {
            respar <- mapply(experiment2parralel,MoreArgs=list(p=p,n=n,pCon=pCon,linear=FALSE, pars = pars),1:numExp)
        }
        
        save.image(paste("./results/experiment2SparseNLP",p,"N",n,".RData", sep=""))
        
        sidlingamDAG[nCounter,pCounter] <- round(mean(respar[7,]),1)
        sidlingamDAGsd[nCounter,pCounter] <- round(sd(respar[7,]),1)
        hdlingamDAG[nCounter,pCounter] <- round(mean(respar[8,]),1)
        hdlingamDAGsd[nCounter,pCounter] <- round(sd(respar[8,]),1)
        hdlingamCPDAG[nCounter,pCounter] <- round(mean(respar[9,]),1)
        hdlingamCPDAGsd[nCounter,pCounter] <- round(sd(respar[9,]),1)
        
        sidicmlDAG[nCounter,pCounter] <- round(mean(respar[4,]),1)
        sidicmlDAGsd[nCounter,pCounter] <- round(sd(respar[4,]),1)
        hdicmlDAG[nCounter,pCounter] <- round(mean(respar[5,]),1)
        hdicmlDAGsd[nCounter,pCounter] <- round(sd(respar[5,]),1)
        hdicmlCPDAG[nCounter,pCounter] <- round(mean(respar[6,]),1)
        hdicmlCPDAGsd[nCounter,pCounter] <- round(sd(respar[6,]),1)
        
        
        sidgdsDAG[nCounter,pCounter] <- round(mean(respar[1,]),1)
        sidgdsDAGsd[nCounter,pCounter] <- round(sd(respar[1,]),1)
        hdgdsDAG[nCounter,pCounter] <- round(mean(respar[2,]),1)
        hdgdsDAGsd[nCounter,pCounter] <- round(sd(respar[2,]),1)
        hdgdsCPDAG[nCounter,pCounter] <- round(mean(respar[3,]),1)
        hdgdsCPDAGsd[nCounter,pCounter] <- round(sd(respar[3,]),1)
        
        sidbfDAG[nCounter,pCounter] <- round(mean(respar[10,]),1)
        sidbfDAGsd[nCounter,pCounter] <- round(sd(respar[10,]),1)
        hdbfDAG[nCounter,pCounter] <- round(mean(respar[11,]),1)
        hdbfDAGsd[nCounter,pCounter] <- round(sd(respar[11,]),1)
        hdbfCPDAG[nCounter,pCounter] <- round(mean(respar[12,]),1)
        hdbfCPDAGsd[nCounter,pCounter] <- round(sd(respar[12,]),1)
        
        sidcpcDAGU[nCounter,pCounter] <- round(mean(respar[17,]),1)
        sidcpcDAGUsd[nCounter,pCounter] <- round(sd(respar[17,]),1)
        sidcpcDAGL[nCounter,pCounter] <- round(mean(respar[18,]),1)
        sidcpcDAGLsd[nCounter,pCounter] <- round(sd(respar[18,]),1)
        hdcpcDAG[nCounter,pCounter] <- round(mean(respar[19,]),1)
        hdcpcDAGsd[nCounter,pCounter] <- round(sd(respar[19,]),1)
        hdcpcCPDAG[nCounter,pCounter] <- round(mean(respar[20,]),1)
        hdcpcCPDAGsd[nCounter,pCounter] <- round(sd(respar[20,]),1)
        
        sidpcDAGU[nCounter,pCounter] <- round(mean(respar[13,]),1)
        sidpcDAGUsd[nCounter,pCounter] <- round(sd(respar[13,]),1)
        sidpcDAGL[nCounter,pCounter] <- round(mean(respar[14,]),1)
        sidpcDAGLsd[nCounter,pCounter] <- round(sd(respar[14,]),1)
        hdpcDAG[nCounter,pCounter] <- round(mean(respar[15,]),1)
        hdpcDAGsd[nCounter,pCounter] <- round(sd(respar[15,]),1)
        hdpcCPDAG[nCounter,pCounter] <- round(mean(respar[16,]),1)
        hdpcCPDAGsd[nCounter,pCounter] <- round(sd(respar[16,]),1)
        
        sidgesDAGU[nCounter,pCounter] <- round(mean(respar[21,]),1)
        sidgesDAGUsd[nCounter,pCounter] <- round(sd(respar[21,]),1)
        sidgesDAGL[nCounter,pCounter] <- round(mean(respar[22,]),1)
        sidgesDAGLsd[nCounter,pCounter] <- round(sd(respar[22,]),1)
        hdgesDAG[nCounter,pCounter] <- round(mean(respar[23,]),1)
        hdgesDAGsd[nCounter,pCounter] <- round(sd(respar[23,]),1)
        hdgesCPDAG[nCounter,pCounter] <- round(mean(respar[24,]),1)
        hdgesCPDAGsd[nCounter,pCounter] <- round(sd(respar[24,]),1)
        
        sidrandDAG[nCounter,pCounter] <- round(mean(respar[25,]),1)
        sidrandDAGsd[nCounter,pCounter] <- round(sd(respar[25,]),1)
        hdrandDAG[nCounter,pCounter] <- round(mean(respar[26,]),1)
        hdrandDAGsd[nCounter,pCounter] <- round(sd(respar[26,]),1)
        hdrandCPDAG[nCounter,pCounter] <- round(mean(respar[27,]),1)
        hdrandCPDAGsd[nCounter,pCounter] <- round(sd(respar[27,]),1)
        
    }
}

cat("======== \n")
cat("SID results \n")
cat("======== \n")
cat("GDS\n")
cat(sidgdsDAG, " +- ", sidgdsDAGsd, "\n \n")
cat("BF\n")
cat(sidbfDAG, " +- ", sidbfDAGsd, "\n \n")
cat("ICML\n")
cat(sidicmlDAG, " +- ", sidicmlDAGsd, "\n \n")
cat("LINGAM\n")
cat(sidlingamDAG, " +- ", sidlingamDAGsd, "\n \n")
cat("cons PC\n")
cat("LB:", sidcpcDAGL, " +- ", sidcpcDAGLsd, "\n")
cat("UB:", sidcpcDAGU, " +- ", sidcpcDAGUsd, "\n \n")
cat("PC\n")
cat("LB:", sidpcDAGL, " +- ", sidpcDAGLsd, "\n")
cat("UB:", sidpcDAGU, " +- ", sidpcDAGUsd, "\n \n")
cat("GES\n")
cat("LB:", sidgesDAGL, " +- ", sidgesDAGLsd, "\n")
cat("UB:", sidgesDAGU, " +- ", sidgesDAGUsd, "\n \n")
cat("random\n")
cat(sidrandDAG, " +- ", sidrandDAGsd, "\n\n")


cat("======== \n")
cat("Hamming Dist results DAG \n")
cat("======== \n")
cat("GDS\n")
cat(hdgdsDAG, " +- ", hdgdsDAGsd, "\n \n")
cat("BF\n")
cat(hdbfDAG, " +- ", hdbfDAGsd, "\n \n")
cat("ICML\n")
cat(hdicmlDAG, " +- ", hdicmlDAGsd, "\n \n")
cat("LINGAM\n")
cat(hdlingamDAG, " +- ", hdlingamDAGsd, "\n \n")
cat("CPC\n")
cat(hdcpcDAG, " +- ", hdcpcDAGsd, "\n \n")
cat("PC\n")
cat(hdpcDAG, " +- ", hdpcDAGsd, "\n \n")
cat("GES\n")
cat(hdgesDAG, " +- ", hdgesDAGsd, "\n \n")
cat("random\n")
cat(hdrandDAG, " +- ", hdrandDAGsd, "\n\n")


cat("======== \n")
cat("Hamming Dist results CPDAG \n")
cat("======== \n")
cat("GDS\n")
cat(hdgdsCPDAG, " +- ", hdgdsCPDAGsd, "\n \n")
cat("BF\n")
cat(hdbfCPDAG, " +- ", hdbfCPDAGsd, "\n \n")
cat("ICML\n")
cat(hdicmlCPDAG, " +- ", hdicmlCPDAGsd, "\n \n")
cat("LINGAM\n")
cat(hdlingamCPDAG, " +- ", hdlingamCPDAGsd, "\n \n")
cat("CPC\n")
cat(hdcpcCPDAG, " +- ", hdcpcCPDAGsd, "\n \n")
cat("PC\n")
cat(hdpcCPDAG, " +- ", hdpcCPDAGsd, "\n \n")
cat("GES\n")
cat(hdgesCPDAG, " +- ", hdgesCPDAGsd, "\n \n")
cat("random\n")
cat(hdrandCPDAG, " +- ", hdrandCPDAGsd, "\n\n")


#for latex
#shd
for(nCounter in 1:length(nVec))
{
    for(pCounter in 1:length(pVec))
    {
        cat("$",hdgdsDAG[nCounter,pCounter], " $&$ ", hdbfDAG[nCounter,pCounter] , " $&$ ", hdicmlDAG[nCounter,pCounter], " $&$ ", hdlingamDAG[nCounter,pCounter], " $&$ ",hdpcDAG[nCounter,pCounter], " $&$ ",hdcpcDAG[nCounter,pCounter]," $&$ ",hdgesDAG[nCounter,pCounter], " $&$ ", hdrandDAG[nCounter,pCounter], "$\n", sep = "")
        cat("$",hdgdsCPDAG[nCounter,pCounter], " $&$ ", hdbfCPDAG[nCounter,pCounter] , " $&$ ", hdicmlCPDAG[nCounter,pCounter], " $&$ ", hdlingamCPDAG[nCounter,pCounter], " $&$ ",hdpcCPDAG[nCounter,pCounter], " $&$ ",hdcpcCPDAG[nCounter,pCounter]," $&$ ",hdgesCPDAG[nCounter,pCounter], " $&$ ", hdrandCPDAG[nCounter,pCounter], "$\n", sep = "")
    }
}

#sid
for(nCounter in 1:length(nVec))
{
    for(pCounter in 1:length(pVec))
    {        
        cat("multirow{2}{*}{$",sidgdsDAG[nCounter,pCounter], "$}&multirow{2}{*}{$", sidbfDAG[nCounter,pCounter] , "$}&multirow{2}{*}{$", sidicmlDAG[nCounter,pCounter], "$}&multirow{2}{*}{$", sidlingamDAG[nCounter,pCounter], "$}&$",sidpcDAGL[nCounter,pCounter], "$&$",sidcpcDAGL[nCounter,pCounter],"$&$",sidgesDAGL[nCounter,pCounter], "$ & multirow{2}{*}{$", sidrandDAG[nCounter,pCounter], "$}\n", sep = "")
        cat("&&&&$",sidpcDAGU[nCounter,pCounter], "$&$",sidcpcDAGU[nCounter,pCounter],"$&$",sidgesDAGU[nCounter,pCounter], "$&\n", sep = "")
    }
}
save.image(paste("./results/experiment2SparseNLP",paste(pVec,collapse=""),"N",paste(nVec,collapse=""),".RData", sep=""))
cat("\n")
