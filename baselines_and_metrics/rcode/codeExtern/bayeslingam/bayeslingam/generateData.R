#Code for regenerating data given the parameters in the MoG case

generateData<-function(DAG,N,components,mixtures=2) {
  nvars<-max(components$node)

  Bup <- dagvector2B( DAG ) #Bup is not lower diagonal!!!!!!!!!!!!!!!

  X<-array(0,c(N,nvars)) 

  for (i in 1:nvars)  {
    node<-DAG[i] #we must go through the nodes in causal order
    #cat('node:');print(node)
    #get the component that is used here!
    for (j in (1:length(components$node)) ) {
      if ( components$node[j] == node && all( components$edges[j,] == Bup[node,] )) 
        icomponent <- j
    }

    #now generate the error term
    X[,node]<-generateMixtureData(N,components$mu[icomponent,1:(3*mixtures)])

    coefindex<-1
    for (j in 1:nvars ) { # then go through all possible parents
      if ( components$edges[icomponent,j] == 1 ) {

        X[,node]<-X[,node]+components$mu[icomponent,3*mixtures+coefindex]*X[,j]
        coefindex<-coefindex+1
      }
    }
  }
  X
}

generateMixtureData<-function(N,mu) { #this more efficient maybe
  mueasy<-parchange(mu) #changing to the easier parameters

  pis<-mueasy[seq(from=1,to=length(mu),by=3)]
  mus<-mueasy[seq(from=2,to=length(mu),by=3)]
  sigmasqs<-mueasy[seq(from=3,to=length(mu),by=3)]

  I<-rmulti(N,pis) #draw from which gaussian
  
  X<-rep(NA,N)

  for (i in 1:length(pis)) {
    w<-I==i
    X[w] <- rnorm(sum(w),mean=mus[i],sd=sqrt(sigmasqs[i]))
  }
  X
}

rmulti<-function(N,prob) {
  R<-runif(N)
  X<-rep(NA,N)
  for (i in seq(from=length(prob),to=1,by=-1)) {
    X[R < sum(prob[1:i])] <- i
  }
  X
}




