structIntervDist2MinAdjSet <- function(trueGraph, estGraph, output = FALSE)
    # Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
    # All rights reserved.  See the file COPYING for license terms. 
{
    estGraph <- as(estGraph, "matrix")
    trueGraph <- as(trueGraph, "matrix")
    G <- trueGraph
    Gp <- estGraph
    PathMatrix <- computePathMatrix(G)
    pMp <- computePathMatrix(Gp)
    incorrectInt <- matrix(0,p,p)
    correctInt <- matrix(0,p,p)
    
    p <- dim(trueGraph)[2]
    for(i in 1:p)
    {
        paG <- which(G[,i] == TRUE)
        for(j in 1:p)
        {
            if(i != j) # test the intervention effect from i to j
            {
                paGp <- computeMinAdjSet(Gp,i,j,pMp)
                #show("--------")
                #show(c(i,j))
                #show(paGp)
                #show("--------")
                #paGp <- which(Gp[,i] == TRUE)
                
                finished <- FALSE
                ijGNull <- FALSE
                ijGpNull <- FALSE
                
                if(PathMatrix[i, j] == 0)
                {
                    ijGNull <- TRUE
                }               
                
                if( (sum(paGp==j)==1) )
                {
                    ijGpNull <- TRUE
                }
                
                
                if(ijGpNull & ijGNull)
                {
                    finished <- TRUE
                    correctInt[i,j] <- 1
                }
                
                if(ijGpNull & !ijGNull)
                {
                    incorrectInt[i,j] <- 1
                    finished <- TRUE
                }
                
                if((length(paGp) == 0) & (!finished))
                {
                    AdjMat2 <- G 
                    AdjMat2[i,] <- rep(FALSE,p)
                    if(dSepAdj(AdjMat2, i, j, c()))
                    {
                        finished <- TRUE
                        correctInt[i,j] <- 1
                    } else
                    {
                        finished <- TRUE
                        incorrectInt[i,j] <- 1
                    }
                } else
                {
                    
                    PathMatrix2 <- computePathMatrix2(G,paGp,PathMatrix)                
                    checkAlldSep <- dSepAdji(G,i,paGp,PathMatrix,PathMatrix2)
                    
                    reachableWOutCausalPath <- checkAlldSep$reachableOnNonCausalPath
                    
                }
                if(!finished && setequal(paG,paGp))
                {
                    finished <- TRUE
                    correctInt[i,j] <- 1
                }
                
                if( !finished )
                {
                    if(PathMatrix[i, j] > 0)
                    {
                        chiCausPath <- which( G[i,] & PathMatrix[,j])
                        if(sum(PathMatrix[chiCausPath,paGp])>0)
                        {
                            incorrectInt[i,j] <- 1
                            finished <- TRUE
                        }
                    }
                    
                    if(!finished)
                    {
                        #check whether all non-causal paths are blocked
                        if(reachableWOutCausalPath[j]==1)
                        {
                            incorrectInt[i,j] <-  1
                        } else
                        {
                            correctInt[i,j] <- 1
                        }
                    }
                    
                }
            }
        }
    }
    
    if(output && p < 11)
    {
        show(incorrectInt)
        show(correctInt)
    }
    ress <- list()
    ress$sid <- sum(incorrectInt)
    ress$incorrectMat <- incorrectInt
    return(ress)
}
