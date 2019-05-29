# extendGreedyResult:
# extend the result of the greedy search by searching more of the close DAGS
#
# SYNTAX:
# results <- extendGreedyResult<-function(R,X,mixtures=2,model='MoG' )
#
#
# INPUT:
# R            - result of the greedy search
# mixtures     - integer (typically 2 or 3): how many gaussian mixture comp
# X            - data, where columns are variables, and rows samples
# model        - 'MoG' or 'GL', this are needed for additional components that
#                might have to be calculated
#
#
# OUTPUT:
# results      - similar type of results list as with the other function
#
#
#

extendGreedyResult<-function(R,X,mixtures=2,model='MoG' ) {
    #this is the thresold, if the DAG prob is under this, its neighbours are no
    #no longer searched
    thresold<-0.005

    X <- normalize(X)


    cat('Extending the most likely dag:\n');
    nvars<-max(R$DAGs)
    R$DAGs<-cdag.to.bdag(R$DAGs)


    extended<-rep(FALSE,length(R$loglike))
    extended[which.max(R$loglike)]<-TRUE
    #extended[which.max(R$prob < thresold)]<-TRUE
    #print(round(100*R$prob))
    #print(round(R$loglike))


    
    while( any(!(extended)) ) {

      #cat('h1\n');


      I<-which(extended == FALSE )

      #cat('h2\n');

      if( length(I) > 1 ) {
        I<-I[1];
      }
      #cat('h3\n');

      if ( R$prob[I] < thresold ) {
        extended[I]<-TRUE
        next;
      }

#      cat('creating neighbours:\n');
      R2<-list()
      R2$DAGs<-neighbourDAGs(R$DAGs[,,I])

      #delete the ones that are already present!  
      #######################################################
      #cat('deleting unnecessary neighbours:\n');

      include<-rep(TRUE,dim(R2$DAGs)[3])
      for ( i in index(1,dim(R2$DAGs)[3]) ) {
        for ( j in index(1,length(R$loglike)) ) {
          if ( equal(R2$DAGs[,,i],R$DAGs[,,j] ) ) {
            include[i]<-FALSE
            break;
          }
        }
      }

      R2$DAGs<-R2$DAGs[,,include]

      #########################################################

      R2$components<-R$components
      R2$loglike<-rep(0,dim(R2$DAGs)[3])




# THIS IS HOW ITS CALLED IN GREEDYBAYESLINGAM!!!

#      cat('calculating scores for extensions:\n');

      R2<-lazylogp(R2,mixtures=mixtures, D=X, model=model)

      R$components<-R2$components #dont forget any new ones calculated

      #calculate the extended dags logs
      #notice now that the extended dag is already in the list
      #also some others might be already on the list

#      cat('calculating probs for extensions:\n');
      R2$prob<-exp(R2$loglike - max(R$loglike) )
      R2$prob<-R2$prob/(sum(R2$prob)+sum(exp(R$loglike - max(R$loglike) ) ) )


      extended[I]<-TRUE
      #add them to R$prob and update its probabilities
      I<- which(R2$prob > thresold)

      if ( length(I) > 0 ) {
        #how to add?
        n1<-dim(R$DAGs)[3]
        n2<-length(I)
        R$DAGs<-c(c(R$DAGs),c(R2$DAGs[,,I]))
        dim(R$DAGs)<-c(nvars,nvars,n1+n2)
        R$loglike<-c(R$loglike, R2$loglike[I])
        extended<-c(extended,rep(FALSE,n2))


        R$prob <- exp(R$loglike - max(R$loglike)) # + priorlog)
        R$prob <- R$prob/sum(R$prob)

      }
      cat('DAGS:',length(extended),'extended ratio:',100*sum(extended)/length(extended),'\n')
        #print(round(100*R$prob))
        #print(round(R$loglike))
  }

  R$DAGs<-bdag.to.cdag(R$DAGs);

  R
}