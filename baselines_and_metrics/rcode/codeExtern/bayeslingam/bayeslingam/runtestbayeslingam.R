runtestbayeslingam<-function(algorithm=c(0,1,2,3,4),dir='plotdata',
                             by=1,from=1,pseudo=FALSE) {
# the numbers of algorithms, 1 for bayeslingam, 2 for deal, 3 for lingam
#   updatedataindex();
#   K<-(dataindex-1)
  #K<-12300
  K<- length(dir(sprintf('%s/',dir),pattern="p[[:digit:]]*.Rdata"))

  Ks<-seq(from=from,to=K,by=by)

  for ( i in 1:length(Ks)) {
        load(sprintf('%s/p%i.Rdata',dir,Ks[i]));#loads parameters

        if ( pseudo ) {
          data<-createPseudoData(parameters);
        } else {
          data<-createData(parameters);
        }
        for ( a in algorithm ) { #test with algorithms in algorithm vector!
          results<-list()

          if ( a==0 ) {
            results2<-bayeslingam(data$X,model='GL')
          } else if ( a==1 ) {
            results2<-bayeslingam(data$X,mixtures=2,model='MoG')
          } else if ( a==2 ) {
            results2<-dealer(data$X)
          } else if ( a == 3 ) {
            results2<-lingamer(data$X)
          } else if ( a == 4 ) {
            results2<-pcer(data$X)
          }

          results$prob<-results2$prob #delete component info
          results$loglike<-results2$loglike
          results$DAGs<-results2$DAGs

          results$version<-20090512

          if ( a==0 ) {
            save(results,file=sprintf('%s/a%i.Rdata',dir,Ks[i]))
          } else if ( a==1 ) {
            save(results,file=sprintf('%s/b%i.Rdata',dir,Ks[i]))
          } else if ( a==2 ) {
            save(results,file=sprintf('%s/g%i.Rdata',dir,Ks[i]))
          } else if ( a == 3 ) {
            save(results,file=sprintf('%s/l%i.Rdata',dir,Ks[i]))
          } else if ( a == 4 ) {
            save(results,file=sprintf('%s/c%i.Rdata',dir,Ks[i]))
          }
        }
  }
}