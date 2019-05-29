indtestMutualHsic <- function(X,alpha=0.05, pars = list())
{
    len <- dim(X)[1]
    p <- dim(X)[2]

    xnorm <- list()
    sigma <- rep(0,p)

    for(i in 1:p)
    {
        xnorm[[i]]<-as.matrix(dist(X[,i],method="euclidean",diag=TRUE,upper=TRUE))
        xnorm[[i]]<-xnorm[[i]]^2
    }
    
    # choose median heuristic for bandwidth
    for(i in 1:p)
    {
        sigma[i] <- sqrt(0.5*median(xnorm[[i]][lower.tri(xnorm[[i]],diag=FALSE)]))
        if(sigma[i] == 0)
        {
            sigma[i] <- sqrt(0.5*mean(xnorm[[i]][lower.tri(xnorm[[i]],diag=FALSE)]))
        }
    }
    
  
    ###
    # Compute kernels
    ###
    K <- list()
    for(i in 1:p)
    {
        K[[i]] <- exp(-xnorm[[i]]/(2*sigma[i]^2))
    }


    ###
    # Compute HSIC
    ###    
    oneVec<- rep(1,len)
    Ktmp <- matrix(1,len,len)
    Ltmp <- matrix(1,len,p)
    Mtmp <- rep(0,p)
    for(i in 1:p)
    {
        Ktmp <- Ktmp * K[[i]]
        Ltmp[,i] <- K[[i]] %*% oneVec
        Mtmp[i] <- sum(K[[i]])
    }
    term1 <- 1/(len^2) * sum(Ktmp)
    term2 <- -2/(len^(p+1)) * sum(apply(Ltmp, 1, prod))
    term3 <- 1/(len^(2*p)) * prod(Mtmp)
    
    HSIC_notsoslow <- term1 + term2 + term3
    # show(HSIC_notsoslow)
    
    return(list(p.value = NaN, statistic = HSIC_notsoslow))
    
    
    ###
    # Compute Gamma Approximation
    ###
    
    #   mux<-(crossprod(colSums(LX))-len)/(len*(len-1))
    #   muy<-(crossprod(colSums(LY))-len)/(len*(len-1))
    
    #   mean_h0<-1/len*(1+mux*muy-mux-muy)
    #   var_h0<-(2*(len-4)*(len-5))/(len*(len-1)*(len-2)*(len-3))*1/((len-1)^2)*sum(diag((t(LX)%*%LXc)%*%(t(LX)%*%LXc)))*1/((len-1)^2)*sum(diag((t(LY)%*%LYc)%*%(t(LY)%*%LYc)))
    
    #   a<-(mean_h0^2)/var_h0
    #   b<-len*var_h0/mean_h0
    
    #   critical_value <- qgamma(1-alpha,shape=a,scale=b)
    #   p_value <- 1-pgamma(len*HSIC,shape=a,scale=b)
    #   Tquan<-c(len*HSIC, critical_value, p_value)
}
