library(reticulate)

args = commandArgs(trailingOnly=TRUE)
# get arguments
new.exp_path <- args[1]   # path to the folder where the results have been saved
new.conda_env <- args[2]
new.conda_path <- args[3]  # can use "auto".

# find .RData file
files <- list.files(path = new.exp_path)
filename <- Filter(function(x) grepl(".RData", x, fixed=TRUE), files)[1]
# load .RData
load(paste(new.exp_path, "/", filename, sep=""))

# hack to make sure loading .RData does not overwrite our data
exp_path <- new.exp_path
conda_env <- new.conda_env
conda_path <- new.conda_path

# reticulate
use_condaenv(condaenv = conda_env, conda = conda_path)

# extract from messy R data structure to a clean python dict
# LINGAM
lingam <- list(sid=respar[7,], shdDAG=respar[8,], shdCPDAG=respar[9,], fdr=respar[43,], time=respar[33,])
# ICML
icml <- list(sid=respar[4,], shdDAG=respar[5,], shdCPDAG=respar[6,], fdr=respar[42,], time=respar[32,])
# GDS
gds <- list(sid=respar[1,], shdDAG=respar[2,], shdCPDAG=respar[3,], fdr=respar[41,], time=respar[31,])
# BF
bf <- list(sid=respar[10,], shdDAG=respar[11,], shdCPDAG=respar[12,], fdr=respar[40,], time=respar[34,])
# CPC
cpc <- list(sidu=respar[17,], sidl=respar[18,], shdDAG=respar[19,], shdCPDAG=respar[20,], fdr=respar[45,], time=respar[36,])
# PC
pc <- list(sidu=respar[13,], sidl=respar[14,], shdDAG=respar[15,], shdCPDAG=respar[16,], fdr=respar[44,], time=respar[35,])
# GES
ges <- list(sidu=respar[21,], sidl=respar[22,], shdDAG=respar[23,], shdCPDAG=respar[24,], fdr=respar[46,], time=respar[37,])
# RANDOM
random <- list(sid=respar[25,], shdDAG=respar[26,], shdCPDAG=respar[27,], fdr=respar[47,], time=respar[38,])
# CAM
cam <- list(sid=respar[28,], shdDAG=respar[29,], shdCPDAG=respar[30,], fdr=respar[48,], time=respar[39,])

# create final dict
dicts <- c(dict(lingam, convert=T), dict(icml, convert=T), dict(gds, convert=T), dict(bf, convert=T),
         dict(cpc, convert=T), dict(pc, convert=T), dict(ges, convert=T), dict(random, convert=T), dict(cam, convert=T))
to_save <- py_dict(c("lingam", "icml", "gds", "bf", "cpc", "pc", "ges", "random", "cam"), dicts, convert=T)

# pickle to python readable format
py_save_object(to_save, paste(exp_path, "/python_results.pkl", sep=""))

