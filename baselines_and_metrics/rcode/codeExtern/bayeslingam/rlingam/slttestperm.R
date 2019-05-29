slttestperm <- function( B )
{
    #Hack of JONAS PETERS 2013
    #
    # slttestperm - tests if we can permute B to strict lower triangularity
    #
    # If we can, then we return the permutation in p, otherwise p=0.
    #
    
    # Dimensionality of the problem
    n <- nrow(B)    
    
    # This will hold the permutation
    p <- c()
    
    # Remaining nodes
    remnodes <- 1:n
    
    # Remaining B, take absolute value now for convenience
    Brem <- abs(B)
    # Select nodes one-by-one
    for(ii in 1:n)
    {
        # Find the row with all zeros
        #therow = find(sum(Brem,2)<1e-12);
        if(length(Brem) > 1)
        {
            rowS <- rowSums(Brem)
        } else
        {
            rowS <- Brem
        }
        therow <- which(rowS < 1e-12)
        
        # If empty, return 0
        #if isempty(therow),
        if(length(therow) == 0)
        {
            p <- 0
            return(p)    
        }
        # If more than one, arbitrarily select the first 
        therow <- therow[1]
        
        # If we made it to the end, then great!
        if(ii==n)
        {
            p <- c(p,remnodes)
            return(p)
        }
        
        # Take out that row and that column
        inds <- which((1:(n-ii+1)) != therow)
        Brem <- Brem[inds,inds]
        ### CHECK!!!!
        
        # Update remaining nodes
        p <- c(p,remnodes[therow])
        remnodes <- remnodes[inds]
        
    }
}