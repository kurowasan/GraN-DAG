#normal priors for MoG model
priorMoG<-list(alpha=list(mu=-1.11,var=0.5,cov=0),
            mu=list(mu=0.0,var=1.0),
            gamma=list(mu=log(0.78),var=1.0),
            b=list(mu=0.0,var=1.0) )


logpriorMoG<-function(p,mixtures=2,prior=priorMoG,derivatives=TRUE) {
  res<-0

  if ( derivatives ) {
    G<-rep(NA,length(p))
    H<-array(0,c(length(p),length(p)))
  }
  m<-mixtures*3#each mixture has 3 parameters

  alphaindex<-seq(from=1,to=m,by=3)
  muindex<-seq(from=2,to=m,by=3)
  gammaindex<-seq(from=3,to=m,by=3)

  if ( length(p) > m ) {
    bindex<-seq(from=m+1,to=length(p),by=1)
  } else {
    bindex<-c()
  }

  alphas<-p[alphaindex]
  mus<-p[muindex]
  gammas<-p[gammaindex]
  bs<-p[bindex]


# 
  #####alpha########
 
  V<-array(0,c(mixtures,mixtures))

   #print(prior)

  diag(V)<-prior$alpha$var

  #here V is just a diagonal so all could be decomposed but also some
  #correlation might be good
  res<-res+(-length(alphas)/2)*log(2*pi)-(1/2)*log(det(V))-
           (1/2)*t(alphas-prior$alpha$mu)%*%solve(V)%*%(alphas-prior$alpha$mu)

  if ( derivatives ) {
    G[alphaindex]<-(-1)*t(alphas-prior$alpha$mu)%*%solve(V)
    H[alphaindex,alphaindex]<-(-1)*solve(V)
  }

  ##########mu################## 

  res<-res+(-length(mus)/2)*log(2*pi)-(length(mus)/2)*log(prior$mu$var)-
          (1/2)*sum((mus-prior$mu$mu)^2)/prior$mu$var

  if ( derivatives ) {
    G[muindex]<-(-(mus-prior$mu$mu)/prior$mu$var)
    H[cbind(muindex,muindex)]<-rep(-1/prior$mu$var,length(mus))
  }

  #####gamma#######

  res<-res+(-length(gammas)/2)*log(2*pi)-
           (length(gammas)/2)*log(prior$gamma$var)-
           (1/2)*sum((gammas-prior$gamma$mu)^2)/prior$gamma$var

  if ( derivatives ) {
    G[gammaindex]<-(-(gammas-prior$gamma$mu)/prior$gamma$var)
    H[cbind(gammaindex,gammaindex)]<-rep(-1/prior$gamma$var,length(gammas))
  }

  #####b#######

  if (!is.null(bindex)) {#only if b is nonempty
    res<-res+(-length(bs)/2)*log(2*pi)-(length(bs)/2)*log(prior$b$var)-
             (1/2)*sum((bs-prior$b$mu)^2)/prior$b$var

    if ( derivatives ) {
      G[bindex]<-(-(bs-prior$b$mu)/prior$b$var)
      H[cbind(bindex,bindex)]<-rep(-1/prior$b$var,length(bs))
    }
  }


  if ( derivatives ) {
    attr(res,"gradient")<-G
    attr(res,"hessian")<-H
  }

  res
}


# Draws 5 random densities from the prior for the MoG density family
drawMoGPrior<-function() {

  # Uses the 2 component gaussian mixture
  mixtures=2;
  
  # Draws five densities and plots them in the same plot
  for (i in 1:5) {
  
    # This actually draws the parameters from the prior
    alpha<-rnorm(mixtures,mean=priorMoG$alpha$mu,sd=sqrt(priorMoG$alpha$var))
    mu<-rnorm(mixtures,mean=priorMoG$mu$mu,sd=sqrt(priorMoG$mu$var))
    gamma<-rnorm(mixtures,mean=priorMoG$gamma$mu,sd=sqrt(priorMoG$gamma$var))

  # Transforms into proportions and variances
    pii<-exp(alpha)/sum(exp(alpha))
    vars<-exp(-gamma)
  
    # Computes the density on this interval
    t<-seq(from=-10,to=10,by=0.01)
    f<-pii[1]*dnorm(t,mean=mu[1],sd=sqrt(vars[1]))+
       pii[2]*dnorm(t,mean=mu[2],sd=sqrt(vars[2]))


  # Plot  
    if ( i==1 ) {
      plot(t,f,ylim=c(0,0.8),xlim=c(-5,5),type='l')
    } else {
      lines(t,f)
    }
  }
}