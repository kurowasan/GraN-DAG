# derivative:
# calculates numerical derivatives for a function
#
# SYNTAX:
# fval<-derivative( p,f,step=1e-5,...)
#
# INPUT:
# p            - the point at which the derivatives are calculated
# f            - the function
# step         - step for the numerical differences
# ...          - additional arguments to f
#
# OUTPUT:
# fval                  - value of f at p
# attr(fval,'gradient') - gradient
# attr(fval,'hessian')  - hessian

derivative<-function( p,f,step=1e-5,...) {

  k<-length(p)

  #3^k experiment design
  P<-array(NA,c(3^k,k))
  for (i in 0:(3^k-1)) {
    P[i+1,]<-dec.to.trin(i,k);
  }
 
  P<-P-1 #changing the values from 0,1,2 to -1,0,1
  P<-P*step #changing to -step,0,step

  y<-rep(NA,3^k)
  for ( i in 1:(3^k) ) {#calculating the values
    y[i]<-value(f(p+P[i,],...))
  }

  Q<-colProduct(P) #quadratic terms

  fit<-lm(y~P+Q)
  #cat('numerical gradient:\n');

  G<-fit$coefficients[2:(k+1)]

  h<-fit$coefficients[(k+2):length(fit$coefficients)]

  H<-array(0,c(k,k));
  H[lower.tri(H,diag=TRUE)]<-h;
  H<-(H+t(H))

  ret<-value(f(p,...))
  attr(ret,'gradient')<-G
  attr(ret,'hessian')<-H

  ret
}

#chages a integer value x to a vector of the trinary system
dec.to.trin <- function( x, ndigits ) {

  b <- array(NA, dim=c(1, ndigits))
  for(i in 1:ndigits) {
    b[, ndigits-i+1] <- (x %% 3)
    x <- (x %/% 3)
  }
  b[1, ]
}