createData<-function(parameters=NULL,dags=NULL,verbal=2) {
  #N must be in the parameters, also nvars
  #dags includes the dags

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
#  data$parameters$theseed<-theseed

#  data$parameters$N<-N
#  data$parameters$nvars<-nvars

#LOGK must be given!!!!!

  if (  parameters$nvars <= 5 ) { 
    if ( is.null(dags)) {
      dags<-alldags(parameters$nvars)
    }
  
    #get the dag #dag index is drawn each time separately -> 
    #it must be different than in the parameters saved in a file!
    #but it should be correct since the seed was used!
    parameters$dagind<-sample(1:nrow(dags),1) 
    parameters$DAG<-dags[parameters$dagind,]
  
    #if ( is.null(parameters$dagind) ) {
    #  parameters$dagind<-sample(1:nrow(dags),1)
    #}
    dag<-dags[parameters$dagind,]

  } else {
    parameters$dagind<-NULL 
    dag<-randomDAG( parameters$nvars );
    parameters$DAG<-dag
  }


  p<- dag[1:parameters$nvars]
  ip <- iperm(p)
  #coefficients
  Nc <- parameters$nvars*(parameters$nvars-1)/2
  coeffs<- dag[parameters$nvars+(1:Nc)]
  parameters$coeffs <- coeffs * runif(Nc,-2.0,2.0)

  #distribution parameters
  parameters$distribution<-array(NA,c(parameters$nvars,5))

  parameters$distribution[,3] <- parameters$logk 
  #here logk can be a vector or a scalar,

  for ( i in 1:parameters$nvars ) {
      #error term parameters
      parameters$distribution[i,1] <- runif(1,-3,3)#mean
      parameters$distribution[i,2] <- runif(1,0.5,3)#standard deviation
      #final parameters of the actual variables
      parameters$distribution[i,4] <- runif(1,-3,3)#mean
      parameters$distribution[i,5] <- runif(1,0.5,3)#standard deviation
  }
  

  #then create data
  B <- matrix(0,parameters$nvars,parameters$nvars)
  B[lower.tri(B,diag=FALSE)] <- parameters$coeffs

  X <- array(0,c(parameters$nvars,parameters$N))

  for (i in 1:(parameters$nvars) ) {
    for (j in 1:i ) { #B is low diagonal!
      X[i,] <- X[i,] + B[i,j]*X[j,]
    }
    k<-exp(parameters$distribution[p[i],3])
    #the error term
    e <- rnorm(parameters$N)
    e <- sign(e)*(abs(e)^k)

    #standardizing mean=0, sd=0
    e <- e-mean(e)
    e <- e/sd(e)
    #print(e)
    #finally setting the mean and the variance
    e <- e*(parameters$distribution[p[i],2])
    e <- e+(parameters$distribution[p[i],1])

    X[i,] <- X[i,] + e

  }
  # permute variable to correct order
  X<-X[ip,]

  #finally set the sds and means according to drawn values,
  # the B matrix is actually no longer valid!
  #but without this we would now in cases where the mean is > 3
  #that there has to be an arrow in to the variable!
  for (i in 1:(parameters$nvars) ) {
    X[i,]<-X[i,]-mean(X[i,])
    X[i,]<-X[i,]/sd(X[i,])
    X[i,]<-X[i,]*parameters$distribution[i,5]
    X[i,]<-X[i,]+parameters$distribution[i,4]
  }

  data$X<-t(X)
  data$parameters<-parameters

  data
}

createDataSet<-function(dir='data',logks,runs,samples,nvars=2) {
  dir.create(dir)
  dags<-alldags(nvars) #create the dags for now

  index<-1
  for (logk in logks) {
    for (run in runs) {
      for (samplesize in samples) {
        parameters<-list()
        parameters$N<-samplesize
        parameters$logk<-logk
        parameters$nvars<-nvars

        data<-createData(parameters,dags=dags) 
        parameters<-data$parameters

        save(parameters,file=sprintf('%s/p%i.Rdata',dir,index))
        #first save parameters
        index<-index+1
      }#samples
    }#runs
  }#logk
}

createDataSetK<-function(dir='data',logks,runs,samples,nvars=2) {
  #for logk plot
  dir.create(dir)
  dags<-alldags(nvars) #create the dags for now

  index<-1
  for (i in 1:nrow(logks)) {
    logk<-logks[i,]
    for (run in runs) {
      for (samplesize in samples) {
        parameters<-list()
        parameters$N<-samplesize
        parameters$logk<-logk
        parameters$nvars<-nvars

        data<-createData(parameters,dags=dags) 
        parameters<-data$parameters

        save(parameters,file=sprintf('%s/p%i.Rdata',dir,index))
        #first save parameters
        index<-index+1
      }#samples
    }#runs
  }#logk
}