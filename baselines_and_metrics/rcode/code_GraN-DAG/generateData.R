library(mgcv)
library(MASS)
library(parallel)
library(pcalg)
library(kernlab)
library(igraph)
library(parallel)
library(gptk)
library(clue)
library(mboost)
library(fastICA)
library(glmnet)
library(CAM)
library(reticulate)
source("../code/util/computeGaussKernel.R", chdir = TRUE)
source("../code/util_DAGs/dag2cpdagAdj.R", chdir = TRUE)
source("../code/util_DAGs/randomDAG.R", chdir = TRUE)
source("../code/util_DAGs/sampleDataFromG.R", chdir = TRUE)
source("../code/util_DAGs/computeCausOrder.R", chdir = TRUE)

# description: Generate DAG, CPDAG and data either from ER or SF graph

# parse command-line arguments
# order: Rscript generateData.R graph_type n p e nb_exp nb_core conda_env conda_path
args = commandArgs(trailingOnly=TRUE)
graph_type <- args[1]   # ER = Erdos-Renyi, SF = Scale-free
n <- as.numeric(args[2])   # number of samples in training set
p <- as.numeric(args[3])  # number of variables in graph
e <- as.numeric(args[4])   # approximate number of edges desired in graph
nb_exp <- as.numeric(args[5])   # number of DAG to generate
nb_core <- as.numeric(args[6])   # number of CPU cores to use
conda_env <- args[7]
conda_path <- args[8]  # can use "auto"

# load conda environment
use_condaenv(condaenv = conda_env, conda = conda_path)
np <- import("numpy")

pCon <- 2 * e / (p**2 - p) # sparse => e edges
funcType <-"GP" # "GAM" or "GP" (function used in GT model)

# create path and file
if(graph_type == 'ER'){
  name.radical = "data"
} else if (graph_type == 'SF'){
  name.radical = "data_sf"
}
data.name <- paste(name.radical, "_p", p, "_e", e, "_n", n, "_", funcType, sep="")
dataset.path <- paste("./data/", data.name, sep="")
dir.create(dataset.path)

for(i in 1:nb_exp)
{
  # Generate data
  if(graph_type == 'ER'){
    trueG <- as.matrix(randomDAG(p,pCon))
  } else if (graph_type == 'SF'){
    trueG <- as.matrix(scaleFreeDAG(p, e))
  }
    
  trueGCPDAG <- dag2cpdagAdj(trueG)
  X <- sampleDataFromG(n,trueG,funcType=funcType, parsFuncType=list(kap=1,sigmax=1,sigmay=1,output=FALSE), noiseType="normalRandomVariances", parsNoise=list(noiseExp=1,varMin=1,varMax=2))
  # parsFuncType explanation
  # B is used only when funcType="linear"
  # kap has an effect only when funcType="GAMGP"
  # sigmax
  
  # parsNoise explanation:
  # noiseVar <- runif(p,parsNoise$varMin,parsNoise$varMax)
  # ran <- rnorm(n)
  # noisetmp <- (0.2*sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)  # if not source node
  # noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)  # if source node
  # (noiseExpVarMin seems to appear only in sampleDataFromGLinear)
  
  # Save DAG, CPDAG and data
  np$save(paste(dataset.path, "/DAG", i, sep=""), trueG)
  np$save(paste(dataset.path, "/CPDAG", i, sep=""), trueGCPDAG)
  np$save(paste(dataset.path, "/data", i, sep=""), X)
}

