# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

#source("../../startups/startupGES.R", chdir = TRUE)
source("../../util_DAGs/randomB.R")
source("../../util_DAGs/randomDAG.R")
source("../../util_DAGs/dag2cpdagAdj.R")
source("../../startups/startupSHD.R", chdir = TRUE)
source("../../startups/startupSID.R", chdir = TRUE)
source("../../startups/startupLINGAM.R", chdir = TRUE)
source("../../startups/startupICML.R", chdir = TRUE)
source("../../startups/startupBF.R", chdir = TRUE)
source("../../startups/startupGDS.R", chdir = TRUE)
source("../../startups/startupPC.R", chdir = TRUE)
source("../../startups/startupScoreSEMIND.R", chdir = TRUE)

stop("The data are not available in this code package.")

load("./Altitude.RData")

resLINGAM <- lingamWrap(cbind(Altitude,Sun,Temp))
cat("LINGAM:\n")
show(resLINGAM$Adj)

resPC <- pcWrap(cbind(Altitude,Sun,Temp),0.01,Inf)
cat("PC:\n")
show(resPC)

resPC <- pcWrap(cbind(Altitude,Sun,Temp),0.01,Inf)
cat("CPC:\n")
show(resPC)

#resGES <- gesWrap(cbind(Altitude,Sun,Temp))
#cat("GES:\n")
#show(resGES$Adj)

# linear
pars <- list(regr.method = train_linear, regr.pars = list(), indtest.method = indtestHsic, indtest.pars = list())
resBF <- BruteForce(cbind(Altitude,Sun,Temp), "SEMIND", pars, output = TRUE)
cat("BF linear:\n")
show(resBF$Adj)

# linear
resICML <- ICML(cbind(Altitude,Sun,Temp),0.05,model=train_linear,indtest=indtestHsic,output= TRUE)
cat("ICML linear:\n")
show(resICML)

# gam
pars <- list(regr.method = train_gam, regr.pars = list(), indtest.method = indtestHsic, indtest.pars = list())
resBFg <- BruteForce(cbind(Altitude,Sun,Temp), "SEMIND", pars, output = TRUE)
cat("BF gam:\n")
show(resBFg$Adj)

# gam
resICML <- ICML(cbind(Altitude,Sun,Temp),0.05,model=train_gam,indtest=indtestHsic,output= TRUE)
cat("ICML gam:\n")
show(resICML)
