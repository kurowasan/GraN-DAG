library(mgcv)
library(pcalg)
library(mboost)
library(CAM)
library(reticulate)

args = commandArgs(trailingOnly=TRUE)
dataset_path <- args[1]   # path to the folder containing the dataset (data1.npy)
dag_path <- args[2]   # path to the folder containing the DAG found by our method (DAG.npy)
exp_path <- args[3]   # path to the folder where the results will be saved
nb_exp <- as.numeric(args[4])   # number of DAG to evaluate
cutOffPVal <- as.numeric(args[5])
conda_env <- args[6]
conda_path <- args[7]  # can use "auto".

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
for (expCounter in 1:nb_exp)
{
    show(paste("Starting exp #", expCounter, sep=""))
    X <- np$load(paste(dataset_path, "/data", expCounter, ".npy", sep=""))  # data generated from the ground truth dag
    # dag <- np$load(paste(dag_path, "/exp", expCounter, "/to-dag/DAG.npy", sep=""))  # dag found by our method
    dag <- np$load(paste(dag_path, "/to-dag/DAG.npy", sep=""))  # dag found by our method
    pruned_dag <- pruning(X, dag, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = cutOffPVal, numBasisFcts = 10), output=TRUE)
    # np$save(paste(exp_path, "/exp", expCounter, "/DAG_pruned.npy" , sep=""), pruned_dag)
    np$save(paste(exp_path, "/pruned_cam/DAG_", cutOffPVal, ".npy" , sep=""), pruned_dag)
}

