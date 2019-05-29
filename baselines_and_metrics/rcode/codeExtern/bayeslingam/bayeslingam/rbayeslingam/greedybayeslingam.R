greedybayeslingam<-function( X,  dags=NULL, verbal=2, mixtures=2,
                         check=FALSE, model='MoG', mcmc=FALSE ) {

  # SYNTAX:
  # r <- bayeslingam( X, dags=NULL, verbal=2, mixtures=2, check =FALSE )
  #
  # INPUT:
  # X        - data matrix, each column is one variable
  # dags     - matrix of possible dags (default: all possible)
  # verbal   - 0=quiet, 1=iter, 2=verbose (default)
  # mixtures - number of components in gaussian mixture (default: 2),
  #            only in MoG model
  # check    - checks log p(replicated data|M) for the best models, 
  #            currently only in MoG case 
  # model    - gaussian mixture ('MoG') (default) or Gaussian Laplace ('GL')
  #
  # OUTPUT:
  # r$DAGs  - list of DAGs, first permutation then connections
  # r$prob  - posterior probabilities of the DAGs
  #
  # NOTE:
  # Needs to have the 'mclust' package installed
  #

  # subtract mean and divide by standard deviation each column separately
  X <- normalize(X)

  # print some statistics
  if ( verbal >= 2 ) {
    print(summary(X))
    print(var(X))
  }

  # let user know we are starting the analysis
  
  # get all DAGs over the variables
  # (current implementation only practical for nvars <= 4)
  N <- dim(X)[1]
  nvars <- dim(X)[2]

  if ( nvars <= 5 ) {
    cat('\n\nConsider exhaustive search running bayeslingam() instead of this greedy search.\n\n')
  }


  if (verbal>=2) cat('Running greedy search bayeslingam:\n')

  cont = TRUE;

  R<-list(components=list(node=c(),edges=array(0,c(0,nvars)) )  );

  #start with an empty BDAG
  BDAG<-array(0,c(nvars,nvars))
  cat(sprintf('(%i) %f %s\n',0,NA,bdag.to.tdag(BDAG)))

  round<-1
  while( cont ) {
    #get the neighbouring DAGs
    R$DAGs<-neighbourDAGs(BDAG)
    R$loglike<-rep(0,dim(R$DAGs)[3])
    #calculate all scores for all dags lazily as 
    # 1. calculating each component only when needed
    # 2. using already calculated components
    R<-lazylogp(R,mixtures=mixtures, D=X, model=model, mcmc=FALSE)

    winning_index=which.max(R$loglike);

    cont= !equal(R$DAGs[,,winning_index],BDAG);
    BDAG<-R$DAGs[,,winning_index];
#     print(BDAG)
#     print(round)
#     print(R$loglike[winning_index])
#     print(bdag.to.tdag(BDAG)[1,1])
# 
#     cat(sprintf('(%i) \n',round) )
#     cat(sprintf('     %f   \n',R$loglike[winning_index]) )
#     cat(sprintf('        %s\n',bdag.to.tdag(BDAG)[1,1] ) )


    cat(sprintf('(%i) %f %s\n',round,R$loglike[winning_index],bdag.to.tdag(BDAG)[1,1] ) )

    round<-round+1
  }

 #cat('R:');print(R)

#   cat('Winning DAG:\n');
#   print(R$DAGs[,,winning_index]);
#   cat('With score:\n');
#   print(R$loglike[winning_index]);
#   cat('With neighbouring DAG scores:\n');
#   print(R$loglike[-winning_index]);

#  priorlog<-rep(1,length(R$loglike))
#  priorlog<-priorlog/sum(priorlog)
#  priorlog<-log(priorlog)

  #set the probs correctly
  R$prob <- exp(R$loglike - max(R$loglike)) # + priorlog)
  R$prob <- R$prob/sum(R$prob)




  R$DAGs<-bdag.to.cdag(R$DAGs);


#  R$prob<-rep(0,length(R$loglike));
#  R$prob[winning_index]<-1;

#   # print out the results...
#   if (verbal>=2) {
# 
#     cat('logp:s for Components:\n')
#     print(cbind(components$node,components$edges,components$score))
# 
#     cat('logp and percentages for DAGS, Mean, and Sd of replicated data:\n')
#     print(cbind(r$DAGs,r$prob,r$loglike,r$repmean,r$repsd))
#   }
# 
#   # return the DAGs and the corresponding probabilities
# 
#   r
  R
}