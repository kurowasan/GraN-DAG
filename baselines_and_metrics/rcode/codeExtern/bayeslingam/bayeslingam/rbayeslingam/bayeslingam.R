#----------------------------------------------------------------------------
# bayeslingam:
# bayesian lingam 
#----------------------------------------------------------------------------

bayeslingam <- function( X,  dags=NULL, verbal=2, mixtures=2,
                         check=FALSE, model='MoG',mcmc=FALSE ) {

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

  if ( nvars==5 || nvars==6 ) {
    cat('\n\nConsider running greedybayeslingam() instead of this exhaustive search.\n\n')
  }
  if ( nvars>=7 ) {
    stop('Too many variables for the exhaustive search. Use greedybayeslingam().\n');
  }


  if (verbal>=2) cat('Creating all dags.\n')

  if (is.null(dags)) dags <- alldags( nvars )
  ndags <- dim(dags)[1]


  if (verbal>=2) cat('Running bayeslingam:\n')


  # let the user know we will be calculating all component scores
  if (verbal>=2) {
    cat('- calculating ',nvars*2^(nvars-1) ,' component logps:\n')  
  }

  # create all components (a component is a node and a parent set)
  components <- allcomponents(nvars)

  # calculate the component scores...
  components <- calculatecomponents( components=components,
                       mixtures=mixtures, D=X, verbal, means=check,
                       model=model, mcmc=mcmc ) 

  # ...and then sum up the component scores to scores for the dags
  if (verbal>=2) cat('\n- sum them to get scores for all ',ndags,' DAGs:\n') 

  loglike <- rep(0,ndags)
  for (i in 1:ndags) {
    loglike[i] <- fastlogp( dags[i,], components )
  }
  # done!
  if (verbal>=1) cat('\n')

  # turn log-likelihood into probability by exponentiating and normalizing.
  # max(loglike) is just for numerical accuracy
  # NOTICE: hear using the uniform prior over the graphs
  # Uniform prior cancels out

  r <- list()

  r$prob <- exp(loglike - max(loglike))
  r$prob <- r$prob/sum(r$prob)

  # this will hold the results
  r$DAGs <- dags
  r$loglike <- loglike
  r$components <- components

  # perform posterior checking
  if ( check  & model == 'MoG') {

    r$repmean <- rep(NA,length(r$prob))
    r$repsd <- rep(NA,length(r$prob))

    if (verbal>=2) {
      cat('Performing posterior check with replicated data...\n')
    }

    bootstraps <- 10
    totalsamples <- 2*N
    Ls <- rep(NA,bootstraps)
    for (i in 1:nrow(dags) ) {
      if ( round(100*r$prob[i]) > 0 ) {
        Xrep <- generateData(r$DAGs[i,],bootstraps*N,r$components)
        L <- logp_vectorized(r$DAGs[i,],r$components,Xrep)
        for ( b in 1:bootstraps ) {
          I <- sample(totalsamples,N)
          Ls[b] <- sum(L$likelihood[I])+L$prior
        }
        r$repmean[i] <- mean(Ls)
        r$repsd[i] <- sd(Ls)
      }
    }
  }



  # print out the results...
  if (verbal>=2) {

    cat('logp:s for Components:\n')
    print(cbind(components$node,components$edges,components$score))

    cat('logp and percentages for DAGS, Mean, and Sd of replicated data:\n')
    print(cbind(r$DAGs,r$prob,r$prob2,r$prob3,r$loglike,r$repmean,r$repsd))
  }

  # return the DAGs and the corresponding probabilities

  r
}
