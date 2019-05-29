#logp calculates the logp-value of a component given parameters and
#data y, parent data X, also derivatives by default


logp<-function(p,y,X,model='MoG',derivatives=TRUE) {
  if ( model == 'MoG') {
    R<-logpMoG(p,y,X,derivatives=derivatives)
  } else if ( model == 'GL' ) {
    R<-logpGL(p,y,X,derivatives=derivatives)
  }

  R
}

#THE FOLLOWING IS JUST FOR POSTERIOR PREDICATIVE CHECK THAT IS NOT PRESENT
#YET ON THIS VERSION

#NOT CLEANED OR COMMENTED, ONLY LINES CUT

logp_vectorized<-function(DAG,components,X,mixtures=2) { 
  #this will calculate the logp faster
  nvars<-ncol(X)
  L<-rep(0,nrow(X))
  P<-0 #this is for prior
  Bup <- dagvector2B( DAG )
  #print(DAG)
  #print(Bup)

  for (i in 1:nvars)  {
    #get the component that is used here!
    for (j in (1:length(components$node)) ) {
      #print(components$node[j])
      #print(Bup[i,])
      #print( components$edges[j,] )
      if ( components$node[j] == i && all( components$edges[j,] == Bup[i,] )) 
        icomponent <- j
    }
    pa<-which( components$edges[icomponent,] != 0 ) 
    Y<-X[,i]
#    cat('Y1:');print(length(Y))
#    cat('pa:');print(pa)
    if (length(pa) >= 1) {
      for (iPa in 1:length(pa)) {
        Y<-Y-components$mu[icomponent,(3*mixtures+iPa)]*X[,pa[iPa]]
      }
    }
    priorv<-logprior(components$mu[icomponent,1:(3*mixtures+length(pa))],
                     mixtures,prior,derivatives=FALSE)

    attributes(priorv)<-NULL
    P<-P+priorv

    Li<-logMixtureGaussian( parchange(components$mu[icomponent,1:(3*mixtures)]),
                            Y,derivatives=FALSE)
    L<-L+Li
  }
  list(likelihood=L,prior=P)
}

logp_parametrized<-function(DAG,components,X,mixtures=2) {
  #this will calculate the logp faster
  nvars<-ncol(X)
  L<-rep(0,nrow(X))
  P<-0 #this is for prior
  Bup <- dagvector2B( DAG )
  #print(DAG)
  #print(Bup)

  for (i in 1:nvars)  {
    #get the component that is used here!
    for (j in (1:length(components$node)) ) {
      #print(components$node[j])
      #print(Bup[i,])
      #print( components$edges[j,] )
      if ( components$node[j] == i && all( components$edges[j,] == Bup[i,] )) 
        icomponent <- j
    }
    pa<-which( components$edges[icomponent,] != 0 ) 
    Y<-X[,i]
#    cat('Y1:');print(length(Y))
#    cat('pa:');print(pa)
    if (length(pa) >= 1) {
      for (iPa in 1:length(pa)) {
        Y<-Y-components$mu[icomponent,(3*mixtures+iPa)]*X[,pa[iPa]]
      }
    }
    priorv<-logprior(components$mu[icomponent,1:(3*mixtures+length(pa))],
                    mixtures,prior)
    #print(priorv)
    attributes(priorv)<-NULL
    P<-P+priorv

    Li<-logMixtureGaussian( parchange(components$mu[icomponent,1:(3*mixtures)]),
                            Y,derivatives=FALSE)
    L<-L+Li
  }
  list(likelihood=L,prior=P)
}