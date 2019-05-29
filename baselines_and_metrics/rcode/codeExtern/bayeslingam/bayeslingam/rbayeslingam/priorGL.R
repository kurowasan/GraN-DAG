#normal priors for the alternative model GL
erf<-function(x) {
  2*pnorm(sqrt(2)*x)-1
}

priorGL<-list( p1=list(mu=0.0,var=1.0),
                 p2=list(mu=0.0,var=1.0),
                 b=list(mu=0.0,var=1.0) )
            #this b is for the coefficients for parents
            #DONT CONFUSE WITH THE OTHER b IN THE MODEL (which is p2 here)


logpriorGL<-function(p,prior=priorGL,derivatives=TRUE) {
  #cat('logprior:')
  res<-0

  if ( derivatives ) {
    G<-rep(NA,length(p))
    H<-array(0,c(length(p),length(p)))
  }

  ###################parameter p1 (aka a)
  res<-res+(-1/2)*log(2*pi)-(1/2)*log(prior$p1$var)-
           (1/2)*(p[1]-prior$p1$mu)^2/prior$p1$var

  if ( derivatives ) {
    G[1]<-(-(p[1]-prior$p1$mu)/prior$p1$var)
    H[1,1]<-(-1)/prior$p1$var
  }

  #########################parameter p2 (aka log(b) aka Beta)

  res<-res+(-1/2)*log(2*pi)-(1/2)*log(prior$p2$var)-
           (1/2)*(p[2]-prior$p2$mu)^2/prior$p2$var

  if ( derivatives ) {
    G[2]<-(-(p[2]-prior$p2$mu)/prior$p2$var)
    H[2,2]<-(-1)/prior$p2$var
  }

  #####b#######
  if ( length(p) > 2 ) {
    bindex<-seq(from=3,to=length(p),by=1)
  } else {
    bindex<-c()
  }
  bs<-p[bindex]

  if ( !is.null(bindex) ) {#only if b is nonempty
    res<-res+(-length(bs)/2)*log(2*pi)-(length(bs)/2)*log(prior$b$var)-
             (1/2)*sum((p[bindex]-prior$b$mu)^2)/prior$b$var

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

# Draws 5 random densities from the prior for the GL density family
drawGLPrior<-function() {

  # Draws five densities and plots them in the same plot
  for (i in 1:5) {

    # This actually draws the parameters from the prior
    p1<-rnorm(1,mean=priorGL$p1$mu,sd=sqrt(priorGL$p1$var))
    p2<-rnorm(1,mean=priorGL$p2$mu,sd=sqrt(priorGL$p2$var))
    p2<-exp(p2)
    p<-c(p1,p2)

    # Computes the normalizing constant
    n<-1
    C<-n*log((-1)*sqrt(pi)*exp((1/4)*p[1]^2/p[2])*
       (-1+erf((1/2)*p[1]/sqrt(p[2])))/sqrt(p[2]) )

    # Computes the density on this interval
    t<-seq(from=-10,to=10,by=0.01)
    sr<-abs(t)
    sr2<-t^2

    Rest<-(-p[1])*sr-p[2]*sr2
    f<-exp(Rest-C)

    # Plot
    if ( i==1 ) {
      plot(t,f,ylim=c(0,1.2),xlim=c(-5,5),type='l')
    } else {
      lines(t,f)
    }
  }
}