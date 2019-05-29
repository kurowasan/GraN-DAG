pcer<-function(X,dags=NULL) {

 # X<-normalize(X)
  nvars<-ncol(X)
#
#  if( is.null(dags) ) {
#    dags<-alldags(nvars)
#  }

  r <- list()
#  r$DAGs <- dags
#  r$loglike<-rep(-Inf,nrow(dags))
#  r$prob <- rep(0,nrow(dags))

  N<-nrow(X)

  significance_level<-0 #these definitions follow the rule of thumb from the README
  if (N >= 1000 ) { #should this be >
    significance_level<-0.01
  } else if ( N >= 100 ) {
    significance_level<-0.05
  } else if ( N >= 50 ) {
    significance_level<-0.1
  } else {
    significance_level<-0.2 #or 0.15?
  }

  #is this necessary
  names(X)<-paste("X",1:ncol(X),sep="") #names for dealer

  Xpc<-make_continuousdata(data.frame(X))


  R<-pc(Xpc,prgt=significance_level)

  Rnothrees<-R
  Rnothrees[R==3]<-0

  while( cycles(Rnothrees) ) {
    cat('Warning, PC-result contained directed cycles.\n');
    cat('Recalculating with a lower significance level:\n');
    significance_level<-significance_level/2
    print(significance_level)
    R<-pc(Xpc,prgt=significance_level)
    Rnothrees<-R
    Rnothrees[R==3]<-0
  }


  print(R)
  #just get the class
  #cat('here\n');print(R)
  apu<-mdag.to.bdag( R )
  #cat('or here\n');print(apu)
  r$DAGs<- bdag.to.cdag( apu )
  #cat('no\n')

  m<-nrow(r$DAGs)
  r$prob<-rep(1/m,m)
  r$loglike<-rep(0,m)
  #now we need to set the probabilities of the dags to the right ones
  #we set equal values to the dags in the equivalence class


#within R matrix
#R[i,j]=0 if there is no edge between the variables
#R[i,j]=2 if there is edge i->j
#R[i,j]=1 if there is edge i<-j
#R[i,j]=3 if there is edge i-j


  #R is now the matrix of coefficients
#   B<- (abs(R$Bprune) > 1e-5)*1

# 
# 
#   for (i in 1:nrow(dags) ) {
#       #cat('i:');print(i)
#      if ( equivalentDAG(dags[i,],R ) ) {
#        r$prob[i]<-1 #set 1 to all equivalent dags
#        r$loglike[i]<-0
#      }
#   }
  #finally adjust the probabilities
#  r$prob<-r$prob/sum(r$prob)
  r
  #R
}

