#logp-value for MoG model, to be maximized by nonlinear optimization

logpMoG<-function(p,y,X,derivatives=TRUE) {
  #X should be a matrix

  mixtures<-(length(p)-ncol(X))/3

  #priorv<-logprior(p=p,mixtures=mixtures,model='MoG',derivatives=derivatives)
  priorv<-logpriorMoG(p=p,mixtures=mixtures,derivatives=derivatives) 


  loglike<-bParSumLogMixtureGaussian(p,y,X,derivatives=derivatives)

  res<-loglike+priorv

  if( derivatives ) {
    attr(res,"gradient")<-attr(loglike,"gradient")+attr(priorv,"gradient")
    attr(res,"hessian")<-attr(loglike,"hessian")+attr(priorv,"hessian")
  } else {
    attr(res,"gradient")<-NULL
    attr(res,"hessian")<-NULL
  }

  res
}