# gives back all the different dags including hidden variable projections (ie.
# every hidden variable has exactly two children):
# if the input is an integer n then the output is all hidden-variable DAGs over
# the empty (no-edges) DAG over the n observed variables
# if a BDAG (over observed variables) is given the output is all hidden-variable
# DAGs over this given BDAG

# need package "gtools" for fct combinations

allhiddenvardags <- function(n=NULL, BDAG=NULL ) {

  O <- BDAG # observed variables
  nvar <- n

  if (is.null(O) & is.null(nvar)) {
    #cat('Input either bdag over observed variables')
    #cat(' or the number of observed variables. \n')
    return( cat('Input either bdag over observed variables or the number of
    observed variables. \n') )
  }
  else if (is.null(O)) O <- array(0,dim=c(nvar,nvar))
  else if (is.null(n)) nvar <- ncol(O)
  else if (ncol(O) != nvar) {
    return( cat('Dimension of bdag-matrix and variables does not agree. \n') )
  }


  library(gtools)
  allpairs <- combinations(n=nvar,r=2)
  npairs <- nrow(allpairs)

  # number of all projections:
  nproj <- 0
  for (i in 1:npairs) {
    nproj <- nproj + choose(npairs,i)
  }

  res <- array(0,dim=c(nvar+npairs,nvar+npairs,nproj))
  cnt <- 1

  for (i in 1:npairs) {
    hv <- combinations(n=npairs,r=i)
    for (j in 1:nrow(hv)) {
      L <- hiddenvar_matrix(allpairs[hv[j,],],nvar)
      res[,,cnt] <- bdaghidd.to.hdag(O,L)
#      res[[cnt]] <- L
      cnt <- cnt + 1
    }
  }

  res
  # res is a 3 dim matrix, each layer in direction of the third dimension is
  # one possible hidden-variable-dag

}


hiddenvar_matrix <- function( hvL, nvar ) {

  # hvL is list of hidden variables, each row containing two nodes which are
  # children of the same hidden var, for example:
  # hvL = 1 2
  #       1 3
  #       2 4
  # means hidden variables between node 1 and 2, between node 1 and 3 and 
  # between node 2 and 4

  hvL <- as.matrix(hvL)
  if (ncol(hvL)==1) hvL <- t(hvL)
  nL <- nrow(hvL)

  L <- array(0,dim=c(nvar,nL))
  for (i in 1:nL) {
    L[hvL[i,],i] <- 1
  }

  L

}