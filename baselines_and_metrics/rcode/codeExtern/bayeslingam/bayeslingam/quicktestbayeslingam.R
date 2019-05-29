# quicktestbayeslingam - some simple quick tests of the method
#
# (Called by 'testbayeslingam')
#

quicktestbayeslingam <- function( nvars=2, N=2000, nongaussiandata=TRUE,
                                  theseed=NULL, drawgraph=FALSE, model='MoG') {


  # Set the seed (random if not specified)
  if (!is.null(theseed) ) {
    set.seed( theseed )
    cat('Using seed: ',theseed,'\n')
  }
  else {
    theseed <- floor(runif(1,0,1)*1000000)
    cat('Using seed: ',theseed,'\n')
    set.seed( theseed )
  }

  
  #----------------- Simulation parameters: ------------------

  # Distributions to use
  if (nongaussiandata) usedistr <- c('unif', 'chisq', 'gengauss')
  else usedistr <- 'gauss'
  
  # Other possible parameters (for the method, e.g.) go here:



  #----------------- Select a random model: ------------------
  
  # Generate all dags over nvars variables
  D <- alldags(nvars)

  # Randomly select one of them
  dagind <- sample(1:nrow(D),1)
  origdag <- D[dagind,]

  # Randomly draw the coefficients in the model
  Nc <- nvars*(nvars-1)/2
  coeffs <- origdag[nvars+(1:Nc)]
  coeffs <- coeffs * runif(Nc,0.5,2.0)

  # Randomly make some coefficients negative? It should be noted though
  # that if we allow negative weights then it is possible for the model
  # to become very close to being unfaithful, and if the data is gaussian
  # then it will be very difficult to identify the correct model equivalence
  # class. Set 'negativecoefficients' to FALSE to avoid this complication.
  negativecoefficients <- FALSE
  if (negativecoefficients) {
    coeffs <- coeffs * sign(rnorm(Nc))
  }
  
  # Randomly select the distributions of the disturbances
  disttype <- rep('',nvars)
  distparams <- list()
  for (i in 1:nvars) {
    disttype[i] <- usedistr[sample(1:length(usedistr),1)]
    distparams[[i]] <- randomparams(disttype[i])
  }

  
  #----------------- Simulate data from the model: ------------------
  
  # Generate data according to this model
  X <- gendata(origdag, nvars, coeffs, disttype, distparams, N)



  #----------------- Run Bayesian Lingam on the data: ---------------  
  
  # Run bayesian lingam
  r <- bayeslingam( t(X), check=FALSE, model=model )


  #-------------------- Displaying the results -------------------
  
  # Now we plot the original model and the top models as predicted
  # by our method. To plot the DAGs we use graphviz.
  
  if ( drawgraph ) {
    plotResults(r)
  }

  cat('\nQuick Results (check up for complete):\n');
  cleanResults( r )

  cat('\nAnd the original was DAG):\n');print(c(origdag))
#  NULL
}


#-----------------------------------------------------------------------------
# gendata:
# generates data according to the DAG linear causal model model
#-----------------------------------------------------------------------------

gendata <- function(origdag, nvars, coeffs, disttype, distparams, N) {

  # create the matrix B (in causal order)
  B <- matrix(0,nvars,nvars)
  B[lower.tri(B,diag=FALSE)] <- coeffs

  # this will hold the data
  X <- matrix(0,nvars,N)

  # this is the variable permutation and inverse permutation
  p <- origdag[1:nvars]
  ip <- iperm(p)
  
  # assign variables in causal order 
  for (i in 1:nvars) {
    for (j in 1:nvars) {
      X[i,] <- X[i,] + B[i,j]*X[j,]
    }
    X[i,] <- X[i,] + sampledata( disttype[p[i]], distparams[[p[i]]], N )
  }

  # permute variable to correct order
  X <- X[ip,]
  X
  
}



#-----------------------------------------------------------------------------
# sampledata:
# randomly samples (simulates) the values of a disturbance variable,
# given its distribution type and parameters
#-----------------------------------------------------------------------------

sampledata <- function( disttype, distparams, N ) {

  m <- distparams[1]
  s <- distparams[2]
  
  if (disttype == 'gauss') {
    x <- rnorm(N, m, s)    
  } else if (disttype == 'unif') {
    x <- runif(N)
    x <- x-mean(x)
    x <- x/sd(x)
    x <- x*s
    x <- x+m
  } else if (disttype == 'chisq') {
    df <- distparams[3]
    x <- rchisq(N,df)
    x <- x-mean(x)
    x <- x/sd(x)
    x <- x*s
    x <- x+m    
  } else if (disttype == 'gengauss') {
    k <- distparams[3]
    x <- rnorm(N)
    x <- sign(x)*(abs(x)^k)
    x <- x-mean(x)
    x <- x/sd(x)
    x <- x*s
    x <- x+m        
  }
  
  x
}
  

#-----------------------------------------------------------------------------
# randomparams:
# select random parameters for the distribution
#-----------------------------------------------------------------------------

randomparams <- function( disttype ) {

  m <- runif(1,-3,3)
  s <- runif(1,0.5,2)
  
  if (disttype == 'gauss') {
    params <- c(m, s)
  } else if (disttype == 'unif') {
    params <- c(m, s)
  } else if (disttype == 'chisq') {
    df <- sample(1:5,1)
    params <- c(m, s, df)
  } else if (disttype == 'gengauss') {
    if (rnorm(1)>0) k <- runif(1,1.3,1.7)
    else k <- runif(1,0.5,0.8)
    params <- c(m, s, k)
  }

  params
  
}



