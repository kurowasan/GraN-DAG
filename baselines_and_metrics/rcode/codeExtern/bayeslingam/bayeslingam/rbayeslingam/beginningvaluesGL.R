# beginningvalues - calculates initial values for numerical optimization
#
# SYNTAX:
# beginningvaluesGL( y, X )
#
# INPUTS:
# y        - the values of the child
# X        - the values of the parents (must be matrix N x number of parents)
#
# OUTPUT:
# p        - vector of starting parameter values
#


beginningvaluesGL<-function(y,X) {


    # the density has two parameters, and the rest are coefficients
    p<-rep(NA,2+ncol(X))
    bindex<-index(3,length(p))

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
  # Alternative density model

    sr<-sum(abs(r)) #the sufficient statistics for the grid search
    sr2<-sum(r^2)
    n<-length(r)

    # GRID APPROACH: JUST CALCULATE THE VALUES ON A GRID
    # the for loop is in case we would like to calculate on a multiple grids
    # this doesnt work very well in any case
    step<-0.2
    p1<-seq(from=-10,to=10,by=step) 
    p2<-seq(from=-10,to=10,by=step)


    Z<-array(-Inf,c(length(p1),length(p2)))

    for ( i in 1:length(p1) ) {
      for ( j in 1:length(p2) ) {
        p[1]<-p1[i]
        p[2]<-p2[j]

        L<-loglikelihood2sufficient(p,sr,sr2,n)+
            logprior(p,model='GL',derivatives=FALSE)

        if (is.finite(L) ) {
          Z[i,j]<-L
        }
      }
    }
    p_i<-which.matrix.max(Z)
    p[1:2]<-c(p1[p_i[1]],p2[p_i[2]])

    #draw a countour of the function
    #later in logp the points of the min routine might be added added to this
    #contour(p1,p2,Z,levels=(max(Z)-seq(from=0,to=10000,by=400)),add=FALSE)
    #wait(5)

  # Return the computed initial values
  p
}

