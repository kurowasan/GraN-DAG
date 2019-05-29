#cleanResults
#Nicely printed results for Bayeslingam results
#Only the best ones printed, with Integer percentages
#Quicktestbayeslingam uses

cleanResults<-function(R) {
  #print(R$DAGs)
  iprob<-round(100*R$prob)

  index<-which( iprob != 0 )
  index<-index[order(iprob[index],decreasing=TRUE)]

  I<-1:nrow(R$DAGs)
  if ( length(index) > 1 ) {
    Ret<-data.frame(cbind(R$DAGs[index,],iprob[index],R$loglike[index]))
  } else {
    Ret<-data.frame(t(c(R$DAGs[index,],iprob[index],R$loglike[index])))
  }

  #names for the columns
  names(Ret)<-c(rep(" ",ncol(R$DAGs)),"iprob","loglike")

  print(Ret)

}