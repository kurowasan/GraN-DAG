sltprune <- function( B )
{
    #Hack of JONAS PETERS 2013
    n <- nrow(B)
    rr <- sort(abs(B), index.return = TRUE)
    
    #[y,ind] = sort(abs(B(:)));
    ind <- rr$ix
    
    for(i in ((n*(n+1)/2):(n*n)))
    {
        
        # Copy original B into Bi
        Bi <- B
        
        # Set 'i' smallest (in absolute value) coefficients to zero
        Bi[ind[1:i]] <- 0
        
        # Try to do permutation
        p <- slttestperm( Bi )
        
        # If we succeeded, then we're done!
        if(any(p != 0))
        {
            Bopt <- B[p,p]
            optperm <- p
            return(list(Bopt = Bopt, optperm = p))
        }
        # ...else we continue, setting one more to zero!
    }
    show("asdsad")
    return(list(Bopt = Bopt, optperm = p))
}