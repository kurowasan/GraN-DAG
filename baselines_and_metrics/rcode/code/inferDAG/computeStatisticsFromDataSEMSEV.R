computeStatisticsFromDataSEMSEV <- function(X,pars)
# Copyright (c) 2010 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    pars$SigmaHat <- cov(X)
    return(pars)
}
