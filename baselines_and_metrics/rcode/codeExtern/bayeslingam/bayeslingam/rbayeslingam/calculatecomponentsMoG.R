# calculatecomponents:
# calculates the bayeslingam scores of a set of components with MoG model
#
# SYNTAX:
# components <- calculatecomponentsMoG( components, mixtures, D, verbal=2,
#                                    means=FALSE, vars=FALSE)
#
# INPUT:
# components   - list of components (such as obtained by allcomponents())
# mixtures     - integer (typically 2 or 3): how many gaussian mixture comp
# D            - data, where columns are variables, and rows samples
# verbal       - how much to print diagnostic information
# means        - output the expectation of the optimal parameters
# vars         - output the variance of the optimal parameters
#
# OUTPUT:
# components   - same as input, but with $scores added, containing log prob

calculatecomponentsMoG <- function( components, mixtures, D, verbal=2,
                                 means=FALSE, vars=FALSE ) {

  # number of variables
  nodes <- ncol(as.matrix(D))

  # number of samples
  N <- nrow(as.matrix(D))

  # initialize array to hold the scores
  components$score <- rep(NA,length(components$node))

  # initialize the arrays holding the means and variances if requested
  if ( means ) {
    components$mu <- array( NA, c(length(components$node),
                                  3*mixtures+max(components$node)-1))
  }
  if ( vars ) {
    components$var <- array(NA,c(length(components$node),
                                 3*mixtures+max(components$node)-1,
                                 3*mixtures+max(components$node)-1))
  }
  
  # calculate the score for each component
  for ( i in 1:length(components$node) ) {

    # this is the data corresponding to the node (the child)
    y <- D[,components$node[i]]

    # index and data of parents
    pa <- which(components$edges[i,] == 1) 
    X <- as.matrix(D[,pa]) 

    # print information on progress?
    if (verbal>=2) {
      cat('[',i,']',' ') 
      cat(printComponent(components$node[i],pa),'\n')
    }
    # set the prior
    # P<-prior


    # calculate the beginning values using Mclust EM-algorithm
    # or by a grid in the GL-case
    x0 <- beginningvaluesMoG( y, X, mixtures )

    # calculate the function value at the beginning point
    # prior for adjusting MoG prior

    #fval <- mf(p=x0,f=logp,y=y,X=X,model='MoG') #
    fval <- mf(p=x0,f=logpMoG,y=y,X=X,derivatives=TRUE) #

    val <- x0

   # pprev<-x0[1:2]
    # The nonlinear minimization algorithm asks for the function value,
    # gradient value, and hessian value in pretty random order, so
    # for speed we try to evaluate each only when strictly necessary

    updateg<-function(p) {
      val<<-p 
      #fval<<-mf(p=p,f=logp,y=y,X=X,model='MoG')#prior=P,,factor=factor
      fval<<-mf(p=p,f=logpMoG,y=y,X=X,derivatives=TRUE)#prior=P,,factor=factor

    }

    #these return the function value, hessian and the gradient
    g<-function(p) {
      updateg(p)
      value(fval)
    }
    hess<-function(p) {
      updateg(p)
      attr(fval,"hessian")
    }
    grad<-function(p) {
      updateg(p)
      attr(fval,"gradient")
    }


    nlminbobject<-nlminb(start=x0,objective=g,gradient=grad,hessian=hess)

    # finally update the function values, in case they are
    # different from what they are currently
    updateg(p=nlminbobject$par)
    p=nlminbobject$par;

    #Convergence check
    if ( nlminbobject$convergence == 0 ) {
      mi<-(-1)*value(fval)
      di<-log(det(attr(fval,"hessian")))

      if ( !is.finite(di) ) {
        if (verbal>=2) cat('Hessian not positive definite!\n')

        #setting the values such that the score will be -Inf
        mi<-(-Inf)
        di<-0
      }
    } else {
      # in this case the hessian eigenvalues might not be positive,
      # and the gradient might not be zero
      if (verbal>=2) cat('Convergence failed!\n');
      if (verbal>=2) print(nlminbobject$message);
      wait(5)
      #print( attr(fval,'gradient') )
      if ( all( eigen(attr(fval,"hessian"),only.values=TRUE)$values > 0 ) &
           all( abs(attr(fval,'gradient')) < 1 ) ) {
        if (verbal>=2) cat('But eigenvalues and gradient quite ok!\n');
        mi<-(-1)*value(fval)
        di<-log(det(attr(fval,"hessian")))

      } else {
        #setting the values such that the score will be -Inf
        mi<-(-Inf)
        di<-0
      }
    }

     if (verbal>=2) {
       cat('Beginning Point for Optimization:\n');print(x0)
       cat('Optimized Point:\n');print(nlminbobject$par)
       cat('Gradient:\n');print(attr(fval,'gradient'))

       cat("Hessian eigenvalues:\n")
       print(eigen(attr(fval,"hessian"),only.values=TRUE)$values)

     }

     # calculate the score (multiply by factorial of mixtures)
     components$score[i] <- (mi-1/2*di+length(x0)/2*log(2*pi))


       #in the mixture gaussian model there are many peaks,
       #so the integral (=score) must be multiplied or in the log
       #form just added
       #anyway this doesn't affect the result, it cancels out
       #inspection necessary only in cases where some have many peaks
       #and other dont
     components$score[i] <- components$score[i] + log(factorial(mixtures))

     if ( means ) {
       #the expected value for the parameters:
       components$mu[i,1:length(minobject$estimate)] <-minobject$estimate
     }

     if ( vars ) {
       #the hessian for the parameters
       components$var[i,1:length(minobject$estimate),
                      1:length(minobject$estimate)] <- solve(minobject$hessian) 
     }
   }

   #return the components list which includes the scores now
   components
}