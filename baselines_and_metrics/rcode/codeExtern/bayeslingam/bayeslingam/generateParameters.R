#for posterior predivative check
#NOT FULLY IMPLEMENTED

generateParameters<-function(DAG,N,components,mixtures=2) {

#the structure has to be 
#1. a dag
#2. the components of that dag
#3. drawn values for each components
  R<-list()
  R$DAG<-DAG

  R<-components #put here only the ones that are needed


  for ( i in length(components$nodes) ) {
    
  }


}

rmultinorm<-function( N, mu, var ) {
  #creates N*length(mu) matrix of values
  V<-sqrtm(var)
  D<-array(mu,c(length(mu),N))
  #print(D)
  E<-array(rnorm(N*length(mu)),c(length(mu),N))
  D<-D+V%*%E

  t(D)
}