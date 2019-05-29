library(mgcv)
library(pcalg)
library(mboost)
library(CAM)
library(reticulate)

args = commandArgs(trailingOnly=TRUE)
dataset_path <- args[1]   # path to the folder containing the dataset (data1.npy)
data_idx <- args[2]
exp_path <- args[3]   # path to the folder where the results will be saved
cutOffPVal <- as.numeric(args[4])
conda_env <- args[5]
conda_path <- args[6]  # can use "auto".

pruning <-
    function(X, G, output = FALSE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10)) 
    {
        p <- dim(G)[1]
        finalG <- matrix(0,p,p)
        for(i in 1:p)
        {
            parents <- which(G[,i]==1)
            lenpa <- length(parents)
            
            if(output)
            {
                cat("pruning variable:", i, "\n")
                cat("considered parents:", parents, "\n")
            }
            
            if(lenpa>0)
            {
                Xtmp <- cbind(X[,parents],X[,i])
                selectedPar <- pruneMethod(Xtmp, k = lenpa + 1, pars = pruneMethodPars, output = output)
                finalParents <- parents[selectedPar]
                finalG[finalParents,i] <- 1
            }
        }
        
        return(finalG)
    }

selGam <-
    function(X,pars = list(cutOffPVal = 0.001, numBasisFcts = 10),output = FALSE,k)
    {
        result <- list()
        p <- dim(as.matrix(X))
        if(p[2] > 1)
        {
            selVec <- rep(FALSE, p[2])
            mod_gam <- CAM:::train_gam(X[,-k],as.matrix(X[,k]),pars)
            pValVec <- summary.gam(mod_gam$model)$s.pv
            if(output)
            {
                cat("vector of p-values:", pValVec, "\n")
            }
            if(length(pValVec) != length(selVec[-k]))
            {
                show("This should never happen (function selGam).")
            }
            selVec[-k] <- (pValVec < pars$cutOffPVal)
        } else
        {
            selVec <- list()
        }
        return(selVec)
    }


use_condaenv(condaenv = conda_env, conda = conda_path)
np <- import("numpy")

# Main
X <- np$load(paste(dataset_path, "/data", data_idx, ".npy", sep=""))  # data generated from the ground truth dag
dag <- np$load(paste(exp_path, "/to-dag/DAG.npy", sep=""))  # dag found by our method
pruned_dag <- pruning(X, dag, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = cutOffPVal, numBasisFcts = 10), output=TRUE)
dir.create(paste(exp_path, "/pruning" , sep=""))
np$save(paste(exp_path, "/pruning/DAG.npy" , sep=""), pruned_dag)

