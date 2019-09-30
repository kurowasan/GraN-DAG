library(CAM)

selGam <-
    function(X,pars = list(cutOffPVal = 0.001, numBasisFcts = 10),output = FALSE,k) {
        result <- list()
        p <- dim(as.matrix(X))
        if(p[2] > 1) {
            selVec <- rep(FALSE, p[2])
            mod_gam <- CAM:::train_gam(X[,-k],as.matrix(X[,k]),pars)
            pValVec <- summary.gam(mod_gam$model)$s.pv
            pValVec <- summary.gam(mod_gam$model)$s.pv
            if(output) {
                cat("vector of p-values:", pValVec, "\n")
            }
            if(length(pValVec) != length(selVec[-k])) {
                show("This should never happen (function selGam).")
            }
            selVec[-k] <- (pValVec < pars$cutOffPVal)
        } else {
            selVec <- list()
        }
        return(selVec)
    }


pruning <-
    function(X, G, output = FALSE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10)) {
        p <- dim(G)[1]
        finalG <- matrix(0,p,p)
        for(i in 1:p) {
            parents <- which(G[,i]==1)
            lenpa <- length(parents)
            
            if(output) {
                cat("pruning variable:", i, "\n")
                cat("considered parents:", parents, "\n")
            }
            
            if(lenpa>0) {
                Xtmp <- cbind(X[,parents],X[,i])
                selectedPar <- pruneMethod(Xtmp, k = lenpa + 1, pars = pruneMethodPars, output = output)
                finalParents <- parents[selectedPar]
                finalG[finalParents,i] <- 1
            }
        }
        
        return(finalG)
    }

dataset <- read.csv(file='{PATH_DATA}', header=FALSE, sep=",")
dag <- read.csv(file='{PATH_DAG}', header=FALSE, sep=",")
set.seed(42)
pruned_dag <- pruning(dataset, dag, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = {CUTOFF}, numBasisFcts = 10), output={VERBOSE})
write.csv(as.matrix(pruned_dag), row.names = FALSE, file = '{PATH_RESULTS}')
