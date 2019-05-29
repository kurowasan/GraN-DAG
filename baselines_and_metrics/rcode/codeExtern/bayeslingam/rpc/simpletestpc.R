# simpletestpc - PC test using completely random parameters
#
# SYNTAX:
# simpletestpc()
#
# What?
# Randomly selects testing parameter values, generates data, and
# runs PC to estimate the graph. Reports and plots the results.
#

simpletestpc <- function( pval ) {

  #----------------------------------------------------------------------
  # 1. GENERATE A MODEL (RANDOMLY SELECT PARAMETERS)
  #----------------------------------------------------------------------

  # Number of variables to use
  dims <- sample(4:6,1)
  cat('Number of dimensions: ',dims,'\n',sep='')

  # Create the network with random weights
  cat('Creating network...')
  B <- matrix(runif(dims*dims),dims,dims)
  B <- B-diag(diag(B))
  B[upper.tri(B)] <- 0
  B[B<0.6] <- 0
  B <- (B>0)*matrix(rnorm(dims*dims),dims,dims)
  disturbancestd <- runif(dims)+1
  c <- 2*rnorm(dims) 
  cat('Done!\n')


  #----------------------------------------------------------------------
  # 2. GENERATE DATA FROM THE MODEL
  #----------------------------------------------------------------------

  cat('Generating data...');

  # Number of data vectors
  samples <- 10000 

  # This generates the disturbance variables, which are mutually 
  # independent, and non-gaussian
  S <- matrix(rnorm(dims*samples),dims,samples)
  
  # This normalizes the disturbance variables to have the 
  # appropriate scales
  S <- S/(sqrt(colMeans(t(S)^2))/disturbancestd)

  # Now we generate the data one component at a time
  Xorig <- matrix(0,dims,samples)
  for (i in 1:dims) {
    Xorig[i,] <- B[i,]%*%Xorig + S[i,] + c[i]
  }
  
  # Select a random permutation because we do not assume that we know
  # the correct ordering of the variables
  p <- sample(1:dims)
  
  # Permute the rows of the data matrix, to give us the observed data
  X <- Xorig[p,]

  # Permute the rows and columns of the original generating matrix B 
  # so that they correspond to the actual data
  Bp <- B[p,p]

  # Permute the generating disturbance stds so that they correspond to
  # the actual data
  disturbancestdp <- disturbancestd[p]

  # Permute the generating constants so that they correspond to
  # the actual data
  cp <- c[p];

  cat('Done!\n')

  #----------------------------------------------------------------------
  # 3. CALL PC TO DO THE ESTIMATION
  #----------------------------------------------------------------------

  # Call PC
  X <- t(X)
  names(X) <- paste("X",1:ncol(X),sep="") 
  Xpc <- make_continuousdata(data.frame(X))
  res <- pc( Xpc, prgt=pval )
  cat('Output of PC:\n')
  print(res)
  cat('\n')
  
  # Process the output 
  Best <- matrix(0,dims,dims)
  for (i in 1:dims) {
    for (j in 1:dims) {

      if (res[j,i]>=2) Best[i,j] <- 1
      
    }
  }
    
  #----------------------------------------------------------------------
  # 4. PLOT THE ESTIMATED GRAPH AGAINST THE TRUE GRAPH
  #----------------------------------------------------------------------

  # For small dimensions, also display the actual connection matrices
  if (dims<=8) {
    print(Bp)
    print(Best)
  }
  
  # Plotting the graphs
  plotgraph(Bp)
  plotgraph(Best)
  

}
