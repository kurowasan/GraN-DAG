#Performs the mcmc analysis for a component, GL model
#return a straight score for it
#This one has it's own prior definitions

mcmcanalysis <- function( X, y ) {

  # Inputs:
  # X one of the following:
  # - NULL (no parents)
  # - a vector (one parent)
  # - an N-by-M matrix (M>1 parents)
  #
  # y is always a vector with the output values
  #
  
  library('adapt')

  # this function computes the logarithm of the posterior density for the
  # supplied parameters, up to the constant given by the partition function
  # (pass NULL 'cvec' and 'X' if no parents!)
  logscore <- function( a, logb, cvec, X, y ) {

    # check if there are no parents
    if (is.null(cvec)) parents <- FALSE
    else parents <- TRUE
    
    # prior over a: (check!)
    amean <- 0
    asd <- 1

    # prior over b: (check!)
    logbmean <- 0
    logbsd <- 1

    # prior over c: (check!)
    if (parents) {
      cmean <- 0
      csd <- 1
    }
      
    # compute the logarithm of the prior density of the given params
    logpa <- dnorm( a, amean, asd, log=TRUE )
    logplogb <- dnorm( logb, logbmean, logbsd, log=TRUE )
    if (parents) logpcvec <- sum(dnorm( cvec, cmean, csd, log=TRUE ))
    
    # compute the log-likelihood
    b <- exp(logb)
    if (parents) {
      if (length(cvec)>1) r <- y - (X %*% cvec) else r <- y - X*cvec
    }
    else {
      r <- y
    }
    logZ <- log(sqrt(pi/b)) + (a^2)/(4*b) + log(2) +
      pnorm((sqrt(2)*a/(2*sqrt(b))),0,1,lower.tail=FALSE,log.p=TRUE)
    logpy <- -a*abs(r) - b*(r^2) - logZ

    # add log-prior and log-likelihood to get the log-score
    if (parents) logsc <- logpa + logplogb + logpcvec + sum(logpy)
    else logsc <- logpa + logplogb + sum(logpy)
    
    # return the log-score
    logsc

  }

  # wrapper for logscore for use with 'adapt' function
  logscorewrapper <- function( inputvec, X, y ) {
    
    a <- inputvec[1]
    logb <- inputvec[2]
    if (length(inputvec)>2) {
      cvec <- inputvec[3:length(inputvec)]
    }
    else {
      cvec <- NULL
    }
    logs <- logscore(a, logb, cvec, X, y)
    logs
    
  }
  
  # this is the proposal (step)
  proposal <- function( a, logb, cvec,stepa=0.01,steplogb=0.01,stepc=0.01 ) {

    # parameters
#    stepa <- 0.01
#    steplogb <- 0.01
#    stepc <- 0.01
    
    # take the step
    newp <- list()
    newp$a <- a + rnorm(1)*stepa
    newp$logb <- logb + rnorm(1)*steplogb
    if (is.null(cvec)) {
      newp$cvec <- NULL
    } else {      
      newp$cvec <- cvec + rnorm(length(cvec))*stepc
    }
      
    # return the new point
    newp
    
  }
  
  # how many variables in X and how many samples?
  if (is.null(X)) {
    nvars <- 0
    N <- length(y)
  } else {
    if (is.null(dim(X))) {
      nvars <- 1
      N <- length(X)
    } else {
      nvars <- dim(X)[2]
      N <- dim(X)[1]
    }
  }
    
  # check that we have same number of samples in Y
  if (length(y) != N) {
    stop('Number of samples in X and y differ (check orientation of X!)');
  }

  #cat('Centering the data...\n');

  # subtract mean from each column of X
  if (nvars>1) {
    X <- t(t(X) - colMeans(X))
  } else if (nvars>0) {
    X <- X - mean(X)
  }
    
  # subtract mean from y
  y <- y - mean(y)



  # compute initial value of coefficients (just linear least-sq regression)
  if (nvars>0) {
    #print(dim(X))
    #cvec <- solve(t(X) %*% X) %*% t(X) %*% y
    #r <- y - X %*% cvec     
    #Antti changed to this!!!
    #cat('lm:');
    fit<-lm(y ~ X - 1)
    #cat('lm done.');

    r<-fit$residual
    cvec<-fit$coef

    #using std errors for the stepsize
    csteps<-0.5*summary(fit)$coef[,2]
    #cat('lm done done.');
  }
  else {
    cvec <- NULL
    r <- y
    csteps<-NULL
  }
    
 # cat('Initialization done...\n');

  #print(csteps)
  # currently initializing a to zero and b so that it matches the
  # variance of the residual, cvec initialized to least-squares coeffs
  a <- 0
  logb <- log(1/(2*(sd(r)^2)))
  logsc <- logscore( a, logb, cvec, X, y )

  # initialize metropolis
  steps <- 1
  Nsteps <- 5000
  histmat <- matrix(0,3+length(cvec),Nsteps)
  histmat[1,1] <- logsc
  histmat[2,1] <- a
  histmat[3,1] <- logb
  if (nvars>0) histmat[4:(dim(histmat)[1]),1] <- cvec
  
  # start metropolis
  rejected=0;
  total=0;
  for (steps in 2:Nsteps) {
    
    # take a random step
    prop <- proposal( a, logb, cvec, stepc=csteps )

    # evaluate
    anew <- prop$a
    logbnew <- prop$logb
    cvecnew <- prop$cvec    
    logscnew <- logscore( anew, logbnew, cvecnew, X, y )

    # accept step? or just keep old point?
    if ((logscnew > logsc) || (runif(1)<exp(logscnew-logsc))) {
      a <- anew
      logb <- logbnew
      cvec <- cvecnew
      logsc <- logscnew
    } else {
      rejected=rejected+1
    }
    
    histmat[1,steps] <- logsc
    histmat[2,steps] <- a
    histmat[3,steps] <- logb
    total=total+1

    if (nvars>0) histmat[4:(dim(histmat)[1]),steps] <- cvec
    
  }

  doplots <- FALSE
  if (doplots) {
  
    # make sure that there is at least Nwin windows open
    Nwin <- 4
    repeat {
      if (length(dev.list())<Nwin) { X11(); next }
      break
    }

    # plot the objective function
    dev.set(which=2)
    plot(1:Nsteps,histmat[1,],'l')
  
    # plot the optimization trajectory
    dev.set(which=3)
    plot(histmat[2,],histmat[3,],'l')

    # plot the optimization trajectory  
    dev.set(which=4)
    if (nvars>0) {
      plot(histmat[3,],histmat[4,],'l')
    } else {
      plot(histmat[2,],histmat[3,],'l')
    }
    
  }
    
  # select the latter half of the chain, parameters only
  points <- histmat[2:(dim(histmat)[1]),(Nsteps/2+1):Nsteps]

  # do a manual fit using a single gaussian, or use mclust?
  manualnormfit <- TRUE
  if (manualnormfit) {
  
    # fit a gaussian mixture model using mclust
    # (for now, just a singal gaussian, manually)
    meanvec <- rowMeans(points)
    covmat <- cov(t(points))
  
    # compute density of gaussian mixture model at each point
    N <- dim(points)[1]
    Xzm <- points-meanvec
    gaussdens <- (1/(((2*pi)^(N/2))*(det(covmat)^(1/2)))) *
      exp(-0.5*colSums(Xzm*(solve(covmat) %*% Xzm)))  

    # look at histogram of log(gaussdens) minus the logscore
    #dev.set(which=5)
    #plot(log(gaussdens),histmat[1,(Nsteps/2+1):Nsteps])
    
    # approximate logarithm of margina likelihood at max
    logmarglikemcmc <- logscore(meanvec[1], meanvec[2], cvec, X, y ) -
      log((1/(((2*pi)^(N/2))*(det(covmat)^(1/2)))))
    cat('logmarglikemcmc = ')
    print(as.vector(logmarglikemcmc[1]),digits=15)  
     cat('rejected = ', rejected, '/' , total, ' = ', round(rejected/total*100), '%\n')
    cat('And the point estimate:\n');
    print(meanvec)

  }
  else {

    x <- points
    xmclust <- Mclust(t(x))
    gaussdens <- dens( modelName = xmclust$modelName, data = t(x),
                        parameters = xmclust$parameters)
    
    # look at histogram of log(gaussdens) minus the logscore
    dev.set(which=5)
    plot(log(gaussdens),histmat[1,(Nsteps/2+1):Nsteps])

    cat('logmarglikemcmc = ')
    #print(logmarglikemcmc,digits=15)  
    cat('(not computed)\n')
    
  }

  # Finally use 'adapt' to approximate the integral (NOT WORKING)
  doadapt <- FALSE
  if (doadapt) {
    ndim <- dim(points)[1]
    lower <- c(-0.2,-0.85,rep(0.9,nvars))
    upper <- c(0.2,-0.55,rep(1.1,nvars))
    logmarglikeadapt <- adapt( ndim, lower, upper, minpts=100, maxpts=NULL,
                              logscorewrapper, eps=0.01, X, y)
    
    cat('logmarglikeadapt = ')
    print(logmarglikeadapt,digits=15)  
  }

  retval <- as.vector(logmarglikemcmc[1])
  retval
  
}
