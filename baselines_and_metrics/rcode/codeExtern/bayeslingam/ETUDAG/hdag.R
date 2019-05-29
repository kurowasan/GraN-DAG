# --- this file contains several functions useful for hdags (= hidden-variable-
# --- DAGs)

# makes DAG-matrix out of DAG over obsered variables and matrix for hidden
# variables

bdaghidd.to.hdag <- function( BDAG_obs, MAT_hidvar ) {

  O <- BDAG_obs
  nobs <- ncol(O)
  L <- MAT_hidvar

  if (is.null(L)) return(O)

  nvar <- nobs + ncol(L) # present variables
   # n° of vars incl. all possible hidden projection vars:
  ntotal <- nobs + choose(nobs,2)
  HDAG <- array(0,dim=c(ntotal,ntotal))
  HDAG[1:nrow(O),1:nvar] <- cbind(O,L)

  HDAG

}

# separates the obsered and hidden variables in an hdag
# HDAG is a dag with hidden variables created f.ex. by the fct bdaghidd.to.hdag

hdag.to.bdaghidd <- function( HDAG ) {

  n <- varnumber_hdag(HDAG)
  nvar <- n[1]+n[2]
  O <- HDAG[1:n[1],1:n[1]]
  #L <- HDAG[(n[1]+1):n[2],(n[1]+1):n[2]]

  temp <- NULL
  for (j in (n[1]+1):nvar) {
    if (all(HDAG[,j]==0)) temp <- c(temp,j)
  }
  L <- as.matrix(HDAG[1:n[1],setdiff((n[1]+1):nvar,temp)])

  list(O=O, L=L)

}

# gives back the number of observed and hidden variables
# HDAG is a dag with hidden variables created by the function hdag

varnumber_hdag <- function( HDAG) {

  nt <- ncol(HDAG)

  nobs <- -1/2 + sqrt(1/4+2*nt)
  nhid <- nt-nobs

  c(nobs,nhid)

}