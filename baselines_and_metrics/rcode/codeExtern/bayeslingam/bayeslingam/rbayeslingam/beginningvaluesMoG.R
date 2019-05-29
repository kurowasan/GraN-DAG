# beginningvalues - calculates initial values for numerical optimization
#
# SYNTAX:
# beginningvaluesMoG( y, X, mixtures, model )
#
# INPUTS:
# y        - the values of the child
# X        - the values of the parents (must be matrix N x number of parents)
# mixtures - the number of components for the 'mixnorm' model
#
# OUTPUT:
# p        - vector of starting parameter values
#


beginningvaluesMoG<-function(y,X,mixtures) {

  # The mixture-of-gaussians model relies on initial values
  # calculated by the mclust R package
  library('mclust')

  # --- Initialize the parameter vector ----------------------------------


  # total number of parameters of the mixtures
  m <- mixtures*3

  # indices of the various types of parameters
  alphaindex<-seq(from=1,to=m,by=3)
  muindex<-seq(from=2,to=m,by=3)
  gammaindex<-seq(from=3,to=m,by=3)

  # indices of the coefficients
  if ( ncol(X) > 0 ) {
    bindex<-seq(from=m+1,to=m+0+ncol(X),by=1)
  } else {
    bindex<-c()
  }

  # generates the parameter vector
  p<-rep(NA,m+ncol(X))


  # --- Fit the linear coefficients, calculate residuals -----------------

  if ( ncol(X) >= 1 ) {

    fit<-lm(y~X-1)
    p[bindex]<-fit$coefficients
    r<-fit$residual
  } else {
    r<-y
  }

  # --- Here's the beef: calculate starting values -----------------------

  # Gaussian mixture model

    # Compute a fit using exactly 'mixtures' mixtures

    mcfit<-Mclust(r,G=mixtures,warn=FALSE)


    p[muindex]<-mcfit$parameters$mean
    p[gammaindex]<-(-1)*log(mcfit$parameters$var$sigmasq)
    if ( mixtures > 1 ) {

      # The proportions output by Mclust are turned into
      # the alphas needed in our method as we use the
      # parametrization: pi = exp(alpha)/sum(exp(alpha)) 

      # We set the mean of the alpha vector appropriately
      # so that it fits with the prior:

      c<-exp(priorMoG$alpha$mu-mean(log(mcfit$parameters$pro)))
      p[alphaindex]<-log(mcfit$parameters$pro)+log(c)

    } else {
    # Only one mixture component, set it to the prior (will
    # be turned into a proportion of one regardless of its
    # value):

      p[alphaindex]<-priorMoG$alpha$mu
    }

  # Return the computed initial values

  p
}


