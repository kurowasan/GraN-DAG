nzdiaghungarian <- function( W )
{
    # Jonas' quick additional hack to make lingam work in high-dimensional problems
    
    n <- nrow(W)
    S <- matrix(1,n,n)/abs(W)
    
    ####
    ###[c,T]=hungarian(S');
    ###
    c <- as.numeric(solve_LSAP(S))
    
    # Permute W to get Wopt
    Pr <- diag(n)
    Pr <- Pr[,c]
    Wopt <- Pr %*% W
    
    
    # Return the optimal permutation as well
    res <- list()
    res$Wopt <- Wopt
    res$rowp <- iperm(c)
    return(res)
    
}
