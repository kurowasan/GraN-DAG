#logp-value for GL model

logpGL<-function(p,y,X,derivatives=TRUE) {
  #this one is to maximized by nlminb!
  #X should be a matrix
    f<-function(p) { #function definition for the derivative

      priorv<-logpriorGL(p=p,derivatives=FALSE)
      loglike<-loglikelihood2(p,y,X,derivatives=FALSE)

      priorv+loglike
    }

    if (derivatives ) {
      res<-derivative(f=f,p=p,step=1e-3) #numerical derivatives!
    } else {
      res<-f(p)
      attr(res,'gradient')<-NULL
      attr(res,'hessian')<-NULL
    }

  res
}
