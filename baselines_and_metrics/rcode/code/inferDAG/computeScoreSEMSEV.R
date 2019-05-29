computeScoreSEMSEV <- function(State, pars, penFactor, output = FALSE)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    SigmaHat <- pars$SigmaHat
    n <- State$n
    p <- State$p
    # sigmaHatSq <- 1/(pn) sum_k sum(State$eachResVar[k] - overallMean)^2
    # but all the means are by construction of the regression zero.
    sigmaHatSq <- (n * State$sumResVar) / (p*n-1)
    # show(State$sumResVar)
    
    I <- diag(rep(1,p))
    B <- State$B
    numPars <- State$numPars
    
    # the following is only the negloglikelihood up to constants
    negloglik <- n * p/2 * log(2*pi*sigmaHatSq) + n/(2*sigmaHatSq) * sum(diag( t(I-B) %*% (I-B) %*% SigmaHat ))
    penalization <- penFactor * log(n)/2 * numPars
    # minimize this BIC score!
    score <- negloglik + penalization
    if(output)
    {
        show("sigmaHatSq")
        show(sigmaHatSq)
        show("negloglik SEMSEV")
        show(negloglik)
        show("pen SEMSEV")
        show(penalization)
        show(score)
    }
    return(score)
}
