computeGaussKernel <- function(x, sigmay, sigmax)
{
    if(is.matrix(x)==FALSE){
        x<-as.matrix(x)}
    n <- dim(x)[1]
    
    xnorm<-as.matrix(dist(x,method="euclidean",diag=TRUE,upper=TRUE))
    xnorm<-xnorm^2
        
    KX <- sigmay * exp(-xnorm/(2*sigmax^2))
        
    return(KX)
}
