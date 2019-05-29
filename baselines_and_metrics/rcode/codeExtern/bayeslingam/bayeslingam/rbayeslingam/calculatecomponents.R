# calculatecomponents:
# calculates the bayeslingam scores of a set of components
#
# SYNTAX:
# components <- calculatecomponents( components, mixtures, D, verbal=2,
#                                    means=FALSE, vars=FALSE, model='MoG')
#
# INPUT:
# components   - list of components (such as obtained by allcomponents())
# mixtures     - integer (typically 2 or 3): how many gaussian mixture comp
# D            - data, where columns are variables, and rows samples
# verbal       - how much to print diagnostic information
# means        - output the expectation of the optimal parameters
# vars         - output the variance of the optimal parameters
# model        - 'MoG' or 'GL'
#
# OUTPUT:
# components   - same as input, but with $scores added, containing log prob
#

calculatecomponents <- function( components, mixtures, D, verbal=2,
                                 means=FALSE, vars=FALSE, model ='MoG', mcmc=FALSE ) {
  if ( model == 'MoG' ) {
    R<-calculatecomponentsMoG( components=components, mixtures=mixtures, D=D, verbal=verbal,
                                 means=means, vars=vars )
  } else if ( model == 'GL' ) {
    if ( mcmc ) {
      R<-calculatecomponentsGLMCMC( components, D, verbal=verbal, means=means, vars=vars ) 
    } else {
      R<-calculatecomponentsGL( components, D, verbal=verbal, means=means, vars=vars ) 
    }
    #not really implemented
  }
  R
}

#prints one component in a nice form
printComponent<-function(node,pa) {
  if ( !is.null(pa) && length(pa) >= 1 ) {
        pastring<-paste(sprintf(':X%i',pa),collapse='')
        substr(pastring,1,1)<-'|'
  } else {
    pastring<-''
  }

  paste('[X',node,pastring,']',sep='')
}


