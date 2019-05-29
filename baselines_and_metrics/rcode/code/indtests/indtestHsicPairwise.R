indtestHsicPairwise <- function(x,y,alpha=0.05, pars = list())
{
    x <- as.matrix(x)
    y <- as.matrix(y)
    px <- dim(x)[2]
    py <- dim(y)[2]
    result <- matrix(0,px*py,3)
    tmp <- matrix(0,px*py,2)
    counter <- 0
    for(i in 1:px)
    {
        for(j in 1:py)
        {
            counter <- counter + 1
            tmp[counter,] <- c(i,j)
            result[counter,1:3] <- indtest_hsic(x[,i], y[,j], alpha, pars)
        }
    }
    best <- which.min(result[,3])
    final <- result[best,]
    final[2] <- indtest_hsic(x[,tmp[best,1]], y[,tmp[best,2]], alpha/(px*py), pars)[2]
    output <- final
}
