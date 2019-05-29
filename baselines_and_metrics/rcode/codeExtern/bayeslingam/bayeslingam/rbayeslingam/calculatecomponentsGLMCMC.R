# calculatecomponentsGL:
# calculates the bayeslingam scores of a set of components using GL model
# and MCMC sampling, the number of samples is hardcoded in mcmcanalysis
#
# SYNTAX:
# components <- calculatecomponentsGLMCMC( components, mixtures, D, verbal=2,
#                                    means=FALSE, vars=FALSE)
#
# INPUT:
# components   - list of components (such as obtained by allcomponents())
# mixtures     - integer (typically 2 or 3): how many gaussian mixture comp
# D            - data, where columns are variables, and rows samples
# verbal       - how much to print diagnostic information
# means        - output the expectation of the optimal parameters (NOT IMPLEMENTED)
# vars         - output the variance of the optimal parameters (NOT IMPLEMENTED)
#
# OUTPUT:
# components   - same as input, but with $scores added, containing log prob
#

calculatecomponentsGLMCMC<-function( components, D, verbal=2,
                                 means=FALSE, vars=FALSE  ) {
  # number of variables
  nodes <- ncol(as.matrix(D))

  # number of samples
  N <- nrow(as.matrix(D))

  # initialize array to hold the scores
  components$score <- rep(NA,length(components$node))

  # initialize the arrays holding the means and variances if requested
  if ( means ) {
    stop('NOT POSSIBLE WITH GL MODEL!\n')
  }
  if ( vars ) {
    stop('NOT POSSIBLE WITH GL MODEL!\n')
  }
  
  # calculate the score for each component
  for ( i in 1:length(components$node) ) {

    # this is the data corresponding to the node (the child)
    y <- D[,components$node[i]]

    # index and data of parents
    pa <- which(components$edges[i,] == 1) 
    X <- as.matrix(D[,pa]) 

    # print information on progress?
    if (verbal>=2) {
      cat('[',i,']',' ') 
      cat(printComponent(components$node[i],pa),'\n')
    }
    R<-mcmcanalysis(X,y)
    print(R)

    components$score[i]<-R
  }

   #return the components list which includes the scores now
   components
}