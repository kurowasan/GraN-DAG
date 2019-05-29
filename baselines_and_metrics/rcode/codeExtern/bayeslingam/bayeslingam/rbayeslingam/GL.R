# Includes the loglikelihood calculation for the GL model


#with approximative the tanh-logcosh approximation is used instead of abs
#but notice that in this case the constant of the pdf is no longer ok
loglikelihood2<-function(p,y,X,print.level=DPL,derivatives=FALSE,approximative=FALSE,factor=10) {
  # cat('loglikelihood2:');
  n<-length(y)
  
  a<-p[1]
  b<-exp(p[2])

  r<-residual(p[index(3,length(p))],y,X,derivatives=FALSE)

  #
  #the alternative, notice that the constant is only approximate here
  if ( approximative ) {
    sr<-sum(1/factor*log(cosh(factor*r)))
  } else {
    sr<-sum(abs(r))
  }
  sr2<-sum(r^2)


  #better formulation for C, the key is the log.p=TRUE of the pnorm function!
  HardPart<-pnorm(0,mean=a/(2*b),sd=1/sqrt(2*b),log.p=TRUE)

  logC<-n*(log(2)+a^2/(4*b)+1/2*log(pi)-1/2*log(b)+HardPart)

  Rest<-(-a)*sr-b*sr2

  res<-Rest-logC


  if ( derivatives ) {
    k<-length(p)
    H<-array(0,c(k,k))
    D<-rep(NA,k)

 
   # expPart<-n*exp(-(1/4)*a^2/b)  
    #these can be used, but the below formulation
    #might be more accurate with certain parameters
    #HarderPart<-(sqrt(pi)*sqrt(p[2])*2*exp(HardPart))
    #expPartperHarderPart2<-expPart/HarderPart

    expPartperHarderPart<-exp (log(n)-(1/4)*a^2/b -1/2*log(pi)-
                               1/2*log(b)-log(2)-HardPart)


    D[1]<- (-1)*(1/2)*n*a/b+expPartperHarderPart-sr

    D[2]<-(1/4)*n*a^2/b+(1/2)*n-(1/2)*a*expPartperHarderPart-b*sr2

    H[1,1]<-(-1/2)*n/b-(1/2)*a/b*expPartperHarderPart+
            1/n*expPartperHarderPart^2

    H[2,2]<-(-1/4)*n*a^2/b-(1/8)*a^3/b*expPartperHarderPart+
            (1/4)*a*expPartperHarderPart+
            (1/4)/n*a^2*expPartperHarderPart^2-b*sr2

    H[2,1]<-H[1,2]<-(1/2)*n*a/b+(1/4)*a^2/b*expPartperHarderPart-
                    (1/2)*expPartperHarderPart-
                    (1/2)/n*p[1]*expPartperHarderPart^2

    if ( k > 2 ) {
      srnax<-apply(as.vector(r)*X,2,sum)
      #s<-sign(as.vector(r))
      #alternative formulation (constant only approximate!)
      if ( approximative ) {
         s<-tanh(factor*as.vector(r))
      } else {
        s<-sign(as.vector(r))
      }

      ssigrx<-apply(s*X,2,sum)

      D[3:k]<-a*ssigrx+2*b*srnax

      H[1,3:k]<-H[3:k,1]<-ssigrx 
      H[2,3:k]<-H[3:k,2]<-2*b*srnax 

      H[3:k,3:k]<-(-2)*b*t(X)%*%X
    }
    attr(res,'gradient')<-D
    attr(res,'hessian')<-H
  }

  res
}

#likelihood2 with sufficient statistics
#used by the grid search for beginning values
loglikelihood2sufficient<-function(p,sr,sr2,n) {

  a<-p[1]
  b<-exp(p[2])

  logC<-n*(log(2)+a^2/(4*b)+1/2*log(pi)-1/2*log(b)+
          pnorm(0,mean=a/(2*b),sd=1/sqrt(2*b),log.p=TRUE))

  Rest<-(-a)*sr-b*sr2

  Rest-logC
}