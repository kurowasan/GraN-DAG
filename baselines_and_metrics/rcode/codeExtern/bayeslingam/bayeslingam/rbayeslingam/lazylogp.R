#Calculates scores while using as many of the precomputed components in R
# and computing the missing ones, and storing them

lazylogp<-function(R, mixtures, D, model, mcmc=FALSE) {
  #X is the data
  n<-dim(R$DAGs)[3]
  nvars<-dim(R$DAGs)[1]

  #R<-list()
  #R$DAGs<-BDAGs
  
  #R$components<-components

  for ( i in 1:n ) {
    for (node in 1:nvars ) {
      component<-list(node=node, edges=R$DAGs[node,,i]  )
      R<-calculatecomponent( component, R, mixtures, D, model,mcmc=mcmc )

      R$loglike[i]<-R$loglike[i] + R$score
    }
  }

  #forget whatever score was stored
  R$score<-NULL

  R
}


#calculatecomponent searches precomputed components for a component
# and calculates possibly the score of a new one

calculatecomponent<-function(component, R, mixtures, D, model,mcmc=FALSE ) {

  found<-FALSE;
  for ( i in index(1,length(R$components$node)) ) {
    #cat('searching index', i , '\n' );
    if ( R$components$node[i] == component$node && 
         all( R$components$edges[i,] == component$edges ) ) {
      #cat('FOUNDDDDD ALREADY CALCULATED COMPONENT!!!!\n')
      R$score<-R$components$score[i];
      found=TRUE;
      break;
    }
  }
  if ( !found ) {
    #calculate it
    dim(component$edges)<-c(1,ncol(D))

    C<-calculatecomponents( component, mixtures=mixtures, D=D, verbal=0,
                                 means=FALSE, vars=FALSE, model = model, mcmc=mcmc );
    cat('.');
    #then add it
    R$components$node<-c(R$components$node,component$node )
    R$components$score<-c(R$components$score,C$score);
    R$score<-C$score;
    #print(R$components$edges)
    #print(component$edges)
    R$components$edges<-rbind(R$components$edges,component$edges)

  }
  R
}