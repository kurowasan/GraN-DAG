# Copyright (c) 2013-2014  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 


source("../../startups/startupLINGAM.R", chdir = TRUE)
source("../../startups/startupICML.R", chdir = TRUE)
source("../../startups/startupBF.R", chdir = TRUE)
source("../../startups/startupGDS.R", chdir = TRUE)
#source("../../startups/startupGES.R", chdir = TRUE)
source("../../startups/startupPC.R", chdir = TRUE)
source("../../startups/startupScoreSEMIND.R", chdir = TRUE)
pars <- list(regr.method = train_linear, regr.pars = list(), indtest.method = indtestHsic, indtest.pars = list())


# generate data
n <- 500
x <- 2*rnorm(n)^3
y <- x + .5*rnorm(n)^3
z <- y + 1.3*rnorm(n)^3
w <- x + z + .5*rnorm(n)^3
X <- cbind(x,y,z,w)
truth <- cbind(c(0,0,0,0),c(1,0,0,0),c(0,1,0,0),c(1,0,1,0))


# run RESIT (also called ICML in the code)
cat("running RESIT...\n")
resICML <- ICML(X, model = train_linear, indtest = indtestHsic, output = FALSE)


# run brute force (works only up to 4 or 5 nodes) 
cat("running Brute Force Algorithm...\n")
resBF <- BruteForce(X, "SEMIND", pars, output = TRUE)$Adj


# run greedy dag search
cat("running GDS...\n")
resGDS <- GDS(X, "SEMIND", pars, check = "checkUntilFirst", output = TRUE, kvec = c(10000), startAt = "emptyGraph")$Adj


# run greedy equivalence search
cat("The GES is not included in the pcalg package yet (Feb 2014). Therefore it does not work on every computer. This may have changed by now.\n") 
# resGES <- gesWrap(X)$Adj


# run LINGAM
cat("running LINGAM...\n")
resLINGAM <- lingamWrap(X)$Adj


#run PC algorithm 
cat("running PC...\n")
resPC <- pcWrap(X,alpha = 0.01, mmax=Inf)

    
#run conservative PC algorithm 
cat("running conservative PC...\n")
resCPC <- cpcWrap(X, 0.01, mmax = Inf)


cat("=======\n")
cat("RESULTS\n")
cat("=======\n")
cat("Truth:\n")
show(truth)
cat("=======\n")
cat("PC:\n")
show(resPC)
cat("=======\n")
cat("conservative PC:\n")
show(resCPC)
cat("=======\n")
cat("LINGAM:\n")
show(resLINGAM)
cat("=======\n")
cat("RESIT:\n")
show(resICML)
cat("=======\n")
cat("greedy DAG search:\n")
show(resGDS)
cat("=======\n")
cat("brute force :\n")
show(resBF)


