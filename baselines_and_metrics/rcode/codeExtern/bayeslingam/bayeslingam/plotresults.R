#plots the graphs of a bayeslingam result

plotResults<-function(results,origdag=NULL,nodenames=NULL) {

  nvars<-max(results$DAGs[1,])

  plotgraph2start()

  if ( !is.null(origdag) ) {
    plotmydag2( origdag, nvars, 'o', 'original' )
  }

  iprob <- round(results$prob*100)
  #print(iprob)
  s <- sort(results$prob, decreasing=TRUE, index.return=TRUE )
  results$DAGs<-results$DAGs[s$ix,]
  iprob<-iprob[s$ix]
  for (i in (1:(dim(results$DAGs)[1])) ) {
      if ( iprob[i] > 0.5 ) { #draw only dags with nonzero percentage
         if (is.null(nodenames)) {
          plotmydag2( results$DAGs[i,], nvars, sprintf('g%g',i), 
                      sprintf('%g%%',iprob[i] ),nodenames=NULL)
        } else {
           plotmydag2( results$DAGs[i,], nvars, sprintf('g%g',i), 
                       sprintf('%g%%',iprob[i] ), nodenames=nodenames )
        }

      }
  }
  plotgraph2stop()
}
