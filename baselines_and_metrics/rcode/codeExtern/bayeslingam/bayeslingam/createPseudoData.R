#creates pseudodata set

createPseudoData<-function(parameters=NULL,verbal=2) {
  if (!is.null(parameters$theseed) ) {
    set.seed( parameters$theseed )
    cat('Using seed: ',parameters$theseed,'\n')
  }
  else {
    parameters$theseed <- floor(runif(1,0,1)*1000000)
    cat('Using seed: ',parameters$theseed,'\n')
    set.seed( parameters$theseed )
  }

  parameters$version<-20090216

  #parameters$datasetindex<- #in this version this is given!
#  parameters$datasetindex<-sample(20,1) #sample one of the 20 real data sets

  X<-getRealData(parameters$datasetindex,verbal)

  parameters$nvars<-ncol(X)

  if ( verbal >= 2 ) {
    cat('Learning the equivalence class with PC.\n')
  }


  R<-pcer(X)
  if ( verbal >= 2 ) {
    cat('Done.\n')
  }

  parameters$class<-which(R$prob > 1e-5) 
  #these are the indexes of the dags in the equivalence class

  parameters$dagind<-sample(parameters$class,1) 

  parameters$DAG<-R$DAGs[parameters$dagind,]
  #draw one of the dags in the equivalence class

  if (is.null(parameters$N)) {
    parameters$N<-max(c(100,round(nrow(X)/4)))
  }

  Bbin<-cdag.to.bdag(R$DAGs[parameters$dagind,]) 
  #get the B matrix of the chosen dag
  B<-Bbin
  res<-array(NA,c(parameters$N,parameters$nvars))
  
  for (i in 1:nrow(Bbin) ) { 
    #we dont have consider the causal order here, 
    #we are only finding out the parameters
    parents<-which(Bbin[i,] != 0 ) #get the parents of the node
    #cat('parents:\n');print(parents)
    fit<-NULL
    if (length(parents) != 0 ) {
      P<-as.matrix(X[,parents])
      fit<-lm(X[,i]~P) #what about the constant term???
      B[i,parents]<-fit$coefficients[2:(length(parents)+1)]
      #cat('coeffs:');print(P%*%D$parameters$B[i,parents])
      res[,i]<-sample(X[,i]-P%*%B[i,parents],parameters$N) 
      #shuffle the residual and draw the needed samples
    } else {
      res[,i]<-sample(X[,i],parameters$N)
    }
  }
  #NOW WE HAVE THE RESIDUALS FOR ALL VARIABLES,
  # WE ALSO FORGET THE ORIGINAL DATA SET

  #then finally create the data according to B
  #now here we must consider the causal order of the variables

  order<-R$DAGs[parameters$dagind,1:parameters$nvars]
  invorder<-iperm(order)

  Bo<-B[order,order] #changing the order to the causal order!
  res<-t(res[,order])

  X <- array(0,c(parameters$nvars,parameters$N))
  #print(X)
  #print(res)
  #print(B)
  for (i in 1:(parameters$nvars) ) {
    for (j in 1:i ) { #B is low diagonal!
      X[i,] <- X[i,] + Bo[i,j]*X[j,]
    }
    X[i,] <- X[i,] + res[i,]
  }
  # permute variable to correct order
  D<-list()
  D$X<-t(X[invorder,])
  D$parameters<-parameters

  D
}

#this creates the pseudodatasets used for the pseudodataplot

createPseudoDataSets<-function() {
  runs<-100
  for (i in 1:20 ) {
    index<-1
    dir<-sprintf('pseudodata%i',i)
    dir.create(dir)
    for (run in 1:runs) {
      parameters<-list(datasetindex=i)
      data<-createPseudoData(parameters) 
      parameters<-data$parameters
      save(parameters,file=sprintf('%s/p%i.Rdata',dir,index))
      #first save parameters
      index<-index+1
    }#runs
  }
}
