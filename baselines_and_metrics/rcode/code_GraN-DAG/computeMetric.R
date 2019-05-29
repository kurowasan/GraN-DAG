source("../code/util_DAGs/dag2cpdagAdj.R", chdir = TRUE)
source("../code/util_DAGs/sampleDataFromG.R", chdir = TRUE)
source("../code/startups/startupSHD.R", chdir = TRUE)
source("../code/startups/startupSID.R", chdir = TRUE)
library(reticulate)
library(pcalg)

# parse command-line arguments

# order: Rscript computeMetric.R dataset_path data_idx exp_path dag_name conda_env conda_path
args = commandArgs(trailingOnly=TRUE)
dataset_path <- args[1]   # path to the folder containing the ground truth dataset (DAG, CPDAG and data)
data_idx <- args[2]
exp_path <- args[3]   # path to the folder where the results have been saved
dag_name <- args[4]   # name of the files containing the DAG
conda_env <- args[5]
conda_path <- args[6]  # can use "auto".

if (length(args) < 6) {
  stop("Arguments missing", call.=FALSE)
}

use_condaenv(condaenv = conda_env, conda = conda_path)
np <- import("numpy")

sids <- rep(0, 1)
hammingDAGs <- rep(0, 1)
hammingCPDAGs <- rep(0, 1)
fdrs <- rep(0, 1)


# load ground truth DAG, CPDAG
gt.DAG <- np$load(paste(dataset_path, "/DAG", data_idx,".npy", sep=""))
gt.CPDAG <- np$load(paste(dataset_path, "/CPDAG", data_idx, ".npy", sep=""))

# load resulting DAG
exp.DAG <- np$load(paste(exp_path, "/", dag_name, sep=""))

# compute metrics
sids[1] <- structIntervDist(gt.DAG, exp.DAG)$sid
hammingDAGs[1] <- hammingDistance(gt.DAG, exp.DAG)
hammingCPDAGs[1] <- hammingDistance(gt.CPDAG, dag2cpdagAdj(exp.DAG))
fdrs[1] <- computeFDR(gt.DAG, exp.DAG)

# create final dict for python
to_save <- py_dict(c("sid", "shdDAG", "shdCPDAG", "fdr"), list(sids, hammingDAGs, hammingCPDAGs, fdrs), convert=T)
# pickle to python readable format
py_save_object(to_save, paste(exp_path, "/python_results.pkl", sep=""))

sink(paste(exp_path , "/results.txt", sep=""))

# print stuff
cat("======== \n")
cat("SID results \n")
cat("======== \n")
cat(round(mean(sids),1), " +- ", round(sd(sids),1), "\n \n")

cat("======== \n")
cat("Hamming Dist results DAG \n")
cat("======== \n")
cat(round(mean(hammingDAGs),1), " +- ", round(sd(hammingDAGs),1), "\n \n")

cat("======== \n")
cat("Hamming Dist results CPDAG \n")
cat("======== \n")
cat(round(mean(hammingCPDAGs),1), " +- ", round(sd(hammingCPDAGs),1), "\n \n")

cat("======== \n")
cat("FDR results \n")
cat("======== \n")
cat(round(mean(fdrs),1), " +- ", round(sd(fdrs),1), "\n \n")

sink()
