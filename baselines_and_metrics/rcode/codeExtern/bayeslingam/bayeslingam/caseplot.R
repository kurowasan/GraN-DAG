
caseplot<-function() {
  #Create the data, and print some of the components
  D<-createCaseData(list(theseed=680511,nvars=6,N=5000))


  Rb<-bayeslingam(D$X) #this is the exhaustive result with MoG model
  cleanResults(Rb)

  Rg<-greedybayeslingam(D$X) #this is the greedy result with MoG model
  cleanResults(Rg)

  i1<-which.max(Rg$prob)

  Rb3<-bayeslingam(D$X,model='GL',mcmc=FALSE) #this exhaustive with numerical derivatives and GL
  cleanResults(Rb3)

  Rg2<-greedybayeslingam(D$X,model='GL',mcmc=FALSE) #this is the greedy search with MCMC
  cleanResults(Rg2)

  Rlingam<-lingamer(D$X)
  cleanResults(Rlingam)

  Rpc<-pcer(D$X)
  cleanResults(Rlingam)

    #then draw the figure

  DAGS<-array(0,c(6,6,3))

  DAGS[,,1]<-cdag.to.bdag(D$parameters$DAG)

  DAGS[,,2]<-cdag.to.bdag(Rg$DAGs[i1,])

  DAGS[,,3]<-cdag.to.bdag(D$parameters$DAG)

  plotdags(DAGS,subtitles=c('original','candidate 1','candidate 2'))

}