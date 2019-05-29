#plots the mixture distribution and the histogram of the actual data
# or in the GL case the GL model

plotmarginal<-function(p,x=c(0),min=-4,max=4,model ) {
  #par(mfg=c(2,1))
  D<-seq(from=min,to=max, by=0.01);
  f<-rep(NA,length(D))

  if ( model == 'MoG' | model == 'mixnorm' ) {
    q<-parchange(p)#q is the vector of sane parameters

    #can use already coded functions!
    for (i in 1:length(D)) {
      f[i]<-value(sumLogMixtureGaussian(p=q,x=D[i]))
    }

    if ( !is.null(x) ) {
      hist(x,freq=FALSE,ylim=c(0,max(exp(f))),xlim=c(min,max),
              add=FALSE,breaks=41)
    }

    lines(D,exp(f)) #this the total approximated distribution

    #and these are the individual components
    for (j in 1:(length(q)/3) ) {
      fval<-q[3*(j-1)+1]*dnorm(D,mean=q[3*(j-1)+2],sd=sqrt(q[3*(j-1)+3]))
      lines(D,fval,col=2)
    }

  } else if (model=='GL' ) {
    sr<-abs(D)
    sr2<-D^2

    for (i in 1:length(D)) {
      f[i]<-loglikelihood2sufficient(p,sr=sr[i],sr2=sr2[i],n=1)
      #f[i]<-loglikelihood3(p[1:2],y=D[i],X=array(0,c(1,0)))
    }

    #ylim=c(0,max(exp(f))),
    hist(x,freq=FALSE,xlim=c(min,max),breaks=41)


    lines(D,exp(f))
  }
}