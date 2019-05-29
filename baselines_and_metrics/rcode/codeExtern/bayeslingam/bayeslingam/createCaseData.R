createCaseData<-function(parameters,verbal=2) {
  #N must be in the parameters, also nvars
  #dags includes the dags
  #calculate the indegree from the nvars
  D<-list()
  parameters$indegree<-min(c(floor(parameters$nvars/2),3))

  data<-list()

  if (!is.null(parameters$theseed) ) {
    set.seed( parameters$theseed )
    if ( verbal >= 2 ) {
      cat('Using seed: ',parameters$theseed,'\n')
    }
  }
  else {
    parameters$theseed <- floor(runif(1,0,1)*1000000)
    if ( verbal >= 2 ) {
      cat('Using seed: ',parameters$theseed,'\n')
    }
    set.seed( parameters$theseed )
  }
  parameters$version<-20081222


  #these are straight from the lingam package
  parminmax <- c(0.5,1.5) # [min max] standard deviation owing to parents
  errminmax <- c(0.5,1.5) # [min max] standard deviation owing to disturbance

  # Create the network with random weights but according to chosen parameters
  cat('Creating balanced network...')
  #this function is a part of the Lingam package actually!
  res <- randnetbalanced2( parameters$nvars, parameters$indegree, parminmax, errminmax )

  B <- res$B
  disturbancestd <- res$errstd
  c <- 2*rnorm(parameters$nvars) 
  cat('Done!\n')


  #----------------------------------------------------------------------
  # 2. GENERATE DATA FROM THE MODEL
  #----------------------------------------------------------------------

   cat('Generating data...\n');

  # Nonlinearity exponent, selected to lie in [0.5, 0.8] or [1.2, 2.0].
  # (<1 gives subgaussian, >1 gives supergaussian)
  q <- runif(parameters$nvars)*1.6+0.2;    
  q<-sample( c(rep(1,10),seq(0.2,0.9,0.1),seq(1.1,1.8,0.1)),parameters$nvars,replace=TRUE)

  #ind <- which(q>0.8)           
  #q[ind] <- q[ind]+0.4     


  # This generates the disturbance variables, which are mutually 
  # independent, and non-gaussian
  S <- matrix(rnorm(parameters$nvars*parameters$N),parameters$nvars,parameters$N)
  S <- sign(S)*(abs(S)^q);
  
  # This normalizes the disturbance variables to have the 
  # appropriate scales
  S <- S/(sqrt(colMeans(t(S)^2))/disturbancestd)

  # Now we generate the data one component at a time
  Xorig <- matrix(0,parameters$nvars,parameters$N)
   for (i in 1:parameters$nvars) {
    Xorig[i,] <- B[i,]%*%Xorig + S[i,] + c[i]
  }
  
  # Select a random permutation because we do not assume that we know
  # the correct ordering of the variables
  p <- sample(1:parameters$nvars)
  
  # Permute the rows of the data matrix, to give us the observed data
  D$X <- t(Xorig[p,])

  # Permute the rows and columns of the original generating matrix B 
  # so that they correspond to the actual data
  Bp <- B[p,p]

  cat('This is the generating network coefficients:\n');
  print(Bp)

  parameters$q<-q[p]
  parameters$DAG<-bdag.to.cdag(abs(Bp) > 1e-5 )

  D$parameters<-parameters
#  data$X<-t(X)
#  data$parameters<-parameters

#  data
  cat('done.\n');
  D
}

randnetbalanced2 <- function( dims, maxindegree, parminmax, errminmax ) {

  # Number of samples used to estimate covariance structure
  samples <- 10000 
    
  # First, generate errstd
  errstd <- runif(dims)*(errminmax[2]-errminmax[1]) + errminmax[1]

  # Initializations
  X <- matrix(0,dims,samples)
  B <- matrix(0,dims,dims)

  # Go trough each node in turn:
  for (i in 1:dims) {


    indegree<-sample(index(0,maxindegree),1)
    # If 'indegree' is finite, randomly pick that many parents,
    # else, all previous variables are parents
    if (!is.infinite(indegree)) {
      if (i<=indegree) {
          if (i>1) {
            par <- 1:(i-1)
          } else {
            par <- rep(0,0)
          }
      } else { 
          par <- sample(1:(i-1),indegree)
      }
    } else {
      if (i>1) {
        par <- 1:(i-1)
      } else {
        par <- rep(0,0)
      }
    }
    # If node has parents, do the following
    if (length(par)>0) {
  
  # Randomly pick weights
  w <- rnorm(length(par))
  wfull <- matrix(0,i-1,1)
        wfull[par] <- w

  # Calculate contribution of parents
  X[i,] <- t(wfull) %*% X[1:(i-1),]
    
  # Randomly select a 'parents std' 
  parstd <- runif(1)*(parminmax[2]-parminmax[1]) + parminmax[1]
  
  # Scale w so that the combination of parents has 'parstd' std
  scaling <- parstd/sqrt(mean(X[i,]^2));
  w <- w*scaling

  # Recalculate contribution of parents
  wfull <- matrix(0,i-1,1)
        wfull[par] <- w 
  X[i,] <- t(wfull) %*% X[1:(i-1),]
  
  # Fill in B
  B[i,par] <- t(w)

      }
    else {
    # if node has no parents
  
      # Increase errstd to get it to roughly same variance
      parstd <- runif(1)*(parminmax[2]-parminmax[1]) + parminmax[1]
      errstd[i] <- sqrt(errstd[i]^2 + parstd^2)
  
      # Set data matrix to empty
      X[i,] = matrix(0,1,samples)
  
    }
  
    # Update data matrix
    X[i,] = X[i,] + matrix(rnorm(samples),1,samples)*errstd[i]
    
  }

  # Return result
  res <- list()
  res$B <- B
  res$errstd <- errstd
  res
  
}
