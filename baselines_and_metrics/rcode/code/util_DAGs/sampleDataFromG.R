sampleDataFromG <- function(n,G,funcType="GAM", parsFuncType=list(B=randomB(G),kap=1,sigmax=1,sigmay=1,output=FALSE), noiseType="normalRandomVariances", parsNoise=list(noiseExp=1,varMin=1,varMax=2,noiseExpVarMin=2,noiseExpVarMax=4,bound=rep(1,dim(G)[2])))
    # Copyright (c) 2013-2013  Jonas Peters  [peters@stat.math.ethz.ch]
    #                          Jan Ernest    [ernest@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
    # 
    # 
    # Generates n samples according to structural equation models based on the DAG G 
    # with specified function class and given noise distribution. Uses Gaussian processes to sample the 
    # nonlinear functions.
    # 
    # INPUT:
    #   m           number of samples that should be simulated
#   G           adjacency matrix of the full DAG
#   funcType    parameter to choose between different function types. Default is "GAM" which simulates from 
#               an additive model (i.e. each node is an additive function of its parents). "GP" samples from  
#               a fully non-additive function, wheras "GAMGP" interpolates between GAM and GP with parameter kap.
#               If funcType="linear" then the function samples from a linear SEM.
#   parsFuncType 
#       kap     interpolates between a general additive model and a fully nonparametric model.
#               kap == 1 --> GAM with additive noise, kap == 0 --> fully nonparametric with additive noise
#       sigmax
#       sigmay  parameters used to simulate the Gaussian kernels
#   noiseType   specifies the type of additive noise in the model. Default is "normalRandomVariances" which simulates
#               Gaussian noise with random variances. 
#   parsNoise   list of parameters to modify the noise distribution.
#
# OUTPUT:
#   X           n x p matrix containing the n samples of the specified model.    

{
    if(funcType == "linear")
    {
        sampleDataFromGLinear(n=n,G=G,parsFuncType,noiseType,parsNoise)
        
    } else if(funcType == "GAM") 
    {
        parsFuncType$kap = 1
        sampleDataFromGAMGP(n=n,G=G,parsFuncType,noiseType,parsNoise)
        
    } else if(funcType == "GP")
    {
        parsFuncType$kap = 0
        sampleDataFromGAMGP(n=n,G=G,parsFuncType,noiseType,parsNoise)
        
    } else if(funcType == "GAMGP")
    {
        sampleDataFromGAMGP(n=n,G=G,parsFuncType,noiseType,parsNoise)
    } else if(funcType == "Sigmoid")
    {
        sampleDataFromMonotoneSigmoid(n=n,G=G,parsFuncType,noiseType,parsNoise)
    } else 
    {
        stop('This function type does not exist!')
    }
}



sampleDataFromGAMGP <- function(n, G, parsFuncType, noiseType, parsNoise)
    # INPUTS:   n:  number of samples
    #           G:  adjacency matrix of Graph to simulate from
    #           kap linearly interpolates between GAM and fully nonparametric model. -- kap == 1 --> GAM with additive noise, kap == 0 --> fully nonparametric with additive noise
    #           noiseExp: exponent to model deviations from Gaussian noise -- noiseExp == 1 --> Gaussian noise
    # OUTPUTS:  X:      sampled data
    
{
    p <- dim(G)[2]
    X <- matrix(NA,n,p)
    # determine the causal Order which is needed for sampling
    causOrder <- computeCausOrder(G)
    
    if(parsFuncType$output)
    {
        show(causOrder)
    }
    
    # sample noise variances
    noiseVar <- runif(p,parsNoise$varMin,parsNoise$varMax)
    
    # loop through each node according to the causal order
    for(node in causOrder)
    {
        if(parsFuncType$output)
        {
            cat("generating GP for node ", node, "\r")
        }
        paOfNode <- which(G[,node] == 1)

        # simulation of noise at source nodes
        if(length(paOfNode) ==0)
        {
            if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
            {
                ran <- rnorm(n)
                noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
            } else
            {
                error("This noiseType is not implemented yet.")
            }
            X[,node] <- noisetmp
        } else
        {
            nuPa <- length(paOfNode)
            X[,node] <- rep(0,n)
            
            # If kap>0 there is an additive model component
            if(parsFuncType$kap>0)
            {
                for(pa in paOfNode)
                {
                    
                    # sample parameters for Gaussian process 
                    kernPa <- computeGaussKernel(X[,pa],parsFuncType$sigmay,parsFuncType$sigmax)
                    fpa <- mvrnorm(1,rep(0,n),kernPa)
                    X[,node] <- X[,node] + parsFuncType$kap * fpa
                }
            }
            
            # if kap<1 there is a non-additive model component
            if(parsFuncType$kap<1)
            {
                kernAllPa <- computeGaussKernel(X[,paOfNode],parsFuncType$sigmay,parsFuncType$sigmax)
                fAllPa <- mvrnorm(1,rep(0,n),kernAllPa)
                if(parsFuncType$output & (parsFuncType$kap==0))
                {
                    ### INCLUDE ADEQUATE PLOTTING FUNCTION (MOREDIMENSIONAL PLOTS) ###                
                }
                X[,node] <- X[,node] + (1-parsFuncType$kap)*fAllPa 
            }
            
            # Additive noise
            if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
            {
                ran <- rnorm(n)
                noisetmp <- (0.2*sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
            } else
            {
                error("This noiseType is not implemented yet.")
            }
            X[,node] <- X[,node] + noisetmp       
        }
    }
    
    return(X)
}





sampleDataFromMonotoneSigmoid <- function(n, G, parsFuncType, noiseType, parsNoise)
    # INPUTS:   n:  number of samples
    #           G:  adjacency matrix of Graph to simulate from
    #           
    #           
    # OUTPUTS:  X:      sampled data
    #
    # This function samples from modified (MONOTONE) sigmoid function 
    # c*b*(x+a)/(1+abs(b*(x+a))) where the choice of a,b,c is random.
    
{
    p <- dim(G)[2]
    X <- matrix(NA,n,p)
    # determine the causal Order which is needed for sampling
    causOrder <- computeCausOrder(G)
    
    if(parsFuncType$output)
    {
        show(causOrder)
    }
    
    # sample noise variances
    noiseVar <- runif(p,parsNoise$varMin,parsNoise$varMax)
    
    # loop through each node according to the causal order
    for(node in causOrder)
    {
        if(parsFuncType$output)
        {
            cat("generating GP for node ", node, "\r")
        }
        paOfNode <- which(G[,node] == 1)
        # simulation of noise at source nodes
        if(length(paOfNode) ==0)
        {
            if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
            {
                ran <- rnorm(n)
                noisetmp <- (sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
            } else
            {
                error("This noiseType is not implemented yet.")
            }
            X[,node] <- noisetmp
        } else
        {
            nuPa <- length(paOfNode)
            X[,node] <- rep(0,n)
            
            # If kap>0 there is an additive model component
            if(parsFuncType$kap>0)
            {
                for(pa in paOfNode)
                {
                    a.sig <- runif(n=1, min=-2, max=2)
                    bern <- rbinom(1,1,0.5)
                    b.sig <- bern*runif(n=1, min=0.5, max=2) + (1-bern)*runif(n=1, min=-2, max=-0.5)
                    c.sig <- rexp(n=1,rate=4)+1
                    X[,node] <- X[,node] + c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig)))
                    if(parsFuncType$output)
                    {
                        plot(X[,pa],c.sig*(b.sig*(X[,pa]+a.sig))/(1+abs(b.sig*(X[,pa]+a.sig))))
                    }
                }
            }
            
            # Additive noise
            if(noiseType == "normalRandomVariances" || noiseType == "normalRandomVariancesFixedExp")
            {
                ran <- rnorm(n)
                noisetmp <- (0.2*sqrt(noiseVar[node]) * abs(ran))^(parsNoise$noiseExp) * sign(ran)
            } else
            {
                error("This noiseType is not implemented yet.")
            }
            X[,node] <- X[,node] + noisetmp       
        }
    }
    
    return(X)
}






sampleDataFromGLinear <- function(n,G,parsFuncType,noiseType,parsNoise)
{
    p <- dim(G)[2]
    Id <- diag(p)
    
    if(noiseType == "unif")
    {
        N <- (cbind(parsNoise$bound) %*% rep(1,n)) * matrix(runif(n*p) -0.5, nrow = p, ncol = n)
        X <- solve(Id-parsFuncType$B) %*% N
        samples <- t(X)   
    } 
    if(noiseType == "normalRandomVariancesFixedExp")
    {
        N <- matrix(sign(rnorm(n*p)), nrow = p, ncol = n) * matrix(rep(runif(p,parsNoise$varMin,parsNoise$varMax),n),nrow = p, ncol=n) * matrix(abs(rnorm(n*p)), nrow = p, ncol = n)^(matrix(parsNoise$noiseExp,nrow = p, ncol=n))
        X <- solve(Id-parsFuncType$B) %*% N
        samples <- t(X)           
    }
    if(noiseType == "normalRandomVariancesRandomExp")
    {
# varMin seems to be the minimal stand dev not the min var.
        N <- matrix(sign(rnorm(n*p)), nrow = p, ncol = n) * matrix(rep(runif(p,parsNoise$varMin,parsNoise$varMax),n),nrow = p, ncol=n) * matrix(abs(rnorm(n*p)), nrow = p, ncol = n)^(matrix(rep(runif(p,parsNoise$noiseExpVarMin,parsNoise$noiseExpVarMax),n),nrow = p, ncol=n))
        X <- solve(Id-parsFuncType$B) %*% N
        samples <- t(X)          
    }
    if(noiseType == "normalGivenVariances")
    {
        SigmaN <- diag(parsNoise$noiseVariances)
        SigmaX <- solve(Id-parsFuncType$B) %*% SigmaN %*% solve(t(Id - parsFuncType$B))
        samples <- mvrnorm(n, rep(0,p), SigmaX)
        return(samples)
    }
    if(noiseType == "normalRandomVariances")
    {
        N <- matrix(rnorm(n*p), nrow = p, ncol = n) * matrix(rep(runif(p,parsNoise$varMin,parsNoise$varMax),n),nrow = p, ncol=n)
        X <- solve(Id-parsFuncType$B) %*% N
        samples <- t(X)           
    }
    return(samples)
}
