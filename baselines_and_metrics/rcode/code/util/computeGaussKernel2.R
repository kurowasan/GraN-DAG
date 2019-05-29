computeGaussKernel2 <- function(x, sigmay, sigmax)
{
    if(is.matrix(x)==FALSE){
        x<-as.matrix(x)}
    n <- dim(x)[1]
    
    Sigma <- matrix(0,n,n)
    for (i in 1:n)
    {
        for (j in i:n)
        {
            Sigma[i,j] <- sigmay * exp(-0.5*((x[i,]-x[j,])*sigmax %*% t(x[i,]-x[j,])))
        }
    }
    Sigma <- Sigma + t(Sigma) - diag(diag(Sigma))
    return(Sigma)

}
