# Copyright (c) 2013 - 2013  Jan Ernest [ernest@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 

# simulate causal effects in true graph

### NOTE: If computing a CE on the true DAG G_0, one has to specify Xpa instead of X.

computeCausalEffectNonlinear <- function(M, G, i, j, x = mean(X[,i]), X, Fpa, Noi, method = "TrueGraph", interp.method = "Linear")
{
    G[,i] = 0
    
    if(i == j)
    {
        return(x)
    }
    
    Mpath <- computePathMatrix(G)
    
    if(Mpath[i,j] == FALSE)
    {
        if(sum(G[,j]==1)==0)
        {
            if(method == "TrueGraph")
            {
                # Simulation from the (known) Noise variances
                return(rnorm(M)*Noi[j])
            } else if(method == "EstimatedGraph") 
            {
                # Works for any noise distribution
                return(sample(X[,j],M,replace=TRUE))
            } else if(method == "EstimatedGraphGaussianNoise")
            {
                # Assuming Gaussian noise -> Estimation of variance and simulation
                sdest <- sd(X[,j])
                return(rnorm(M)*sdest)
            } else
            {
                stop('The method you specified is not available!')
            }
            isComputed[rpa] <- TRUE
        } else {
            coo <- computeCausOrder(G)
            cO <- coo[which(G[coo,j]==1)[1]]
            return(computeCausalEffectNonlinear(M,G,cO,j,x = rnorm(M)*Noi[cO],X,Fpa,Noi,method=method,interp.method=interp.method))        
        }
    }
    
    p <- dim(G)[1]
    isComputed <- rep(FALSE, p)
    causEffect <- matrix(0,M,p)
    isComputed[i] <- TRUE
    causEffect[,i] <- x
    TT <- vector("list", p)
    
    children <- which(G[i,] == 1 & Mpath[,j] == TRUE)
    
    while(isComputed[j]==FALSE)
    {
        for(node in children)
        {
            ### Compute causal effects coming from direct parents on the path
            all_parents <- which(G[,node] == 1)
            n_parents <- length(all_parents)
            parents_on_path <- which(G[,node] == 1 & Mpath[i,] == TRUE)
            if(sum(isComputed[parents_on_path] == FALSE)==0)
            {
                ### Compute the value of all parents NOT on the path
                parents_not_on_path <- all_parents[!(all_parents%in%parents_on_path)]
                if(length(parents_not_on_path)>0)
                {
                    root_parents <- find_parent_roots(i,node,Mpath,G) 
                    if(!(length(root_parents)==0))
                    {
                        chi <- NULL
                        for(rpa in root_parents)
                        {
                            if(!isComputed[rpa])
                            {
                                if(method == "TrueGraph")
                                {
                                    # Simulation from the (known) Noise variances
                                    causEffect[,rpa] <- rnorm(M)*Noi[rpa]
                                } else if(method == "EstimatedGraph") 
                                {
                                    # Works for any noise distribution
                                    causEffect[,rpa]  <- sample(X[,rpa],M,replace=TRUE)
                                } else if(method == "EstimatedGraphGaussianNoise")
                                {
                                    # Assuming Gaussian noise -> Estimation of variance and simulation
                                    sdest <- sd(X[,rpa])
                                    causEffect[,rpa] <- rnorm(M)*sdest
                                } else
                                {
                                    stop('The method you specified is not available!')
                                }
                                isComputed[rpa] <- TRUE
                            }
                            
                            tmp <- which((G[rpa,]==1) & (Mpath[i,]==0) & (Mpath[,node]==1))
                            chi <- c(tmp[which(!(tmp %in% chi))], chi)
                        }    
                        
                        while(sum(chi!=node)>0)
                        {
                            # while not all causal effects of the parents not on the path are determined
                            for(nod in chi[chi!=node])
                            {
                                all_pa <- which(G[,nod] == 1)
                                n_pa <- length(all_pa)
                                
                                if(sum(isComputed[all_pa] == FALSE)==0)
                                {
                                    if(!isComputed[nod]){
                                        if(method == "TrueGraph")
                                        {
                                            for(pare in all_pa)
                                            {
                                                causEffect[,nod] <- causEffect[,nod] + approx(X[,pare], Fpa[,pare,nod], xout = causEffect[,pare], method = "linear", rule = 1)$y
                                            }
                                            
                                            causEffect[,nod] <- causEffect[,nod] + rnorm(M)*Noi[nod]
                                            
                                        } else if(method == "EstimatedGraph")
                                        {
                                            #TT <- train_gam_causEffect(X[,all_pa],X[,nod],nod,TT)
                                            TT[[nod]] <- train_gam(X=X[,all_pa],y=X[,nod])$model
                                            datafr <- data.frame(causEffect[,all_pa])
                                            nam <- "var2"
                                            if(n_pa > 1)
                                            {
                                                for(ii in 2:n_pa)
                                                {
                                                    nam <- c(nam,paste("var", ii+1, sep = ""))
                                                }
                                            }
                                            names(datafr) <- nam
                                            causEffect[,nod] <- causEffect[,nod] + predict.gam(object = TT[[nod]], newdata = datafr , type = "link") + sample(TT[[nod]]$residuals,M,replace=TRUE)
                                        } else if(method == "EstimatedGraphGaussianNoise")
                                        {
                                            #TT <- train_gam_causEffect(X[,all_pa],X[,nod],nod,TT)
                                            TT[[nod]] <- train_gam(X=X[,all_pa],y=X[,nod])$model
                                            datafr <- data.frame(causEffect[,all_pa])
                                            nam <- "var2"
                                            if(n_pa > 1)
                                            {
                                                for(ii in 2:n_pa)
                                                {
                                                    nam <- c(nam,paste("var", ii+1, sep = ""))
                                                }
                                            }
                                            names(datafr) <- nam
                                            causEffect[,nod] <- causEffect[,nod] + predict.gam(object = TT[[nod]], newdata = datafr , type = "link") + rnorm(M)*sd(TT[[nod]]$residuals)
                                        } else
                                        {
                                            stop('The method you specified is not available')
                                        }
                                    }
                                    chi <- chi[!chi==nod]
                                    tmp <- which((G[nod,]==1) & (Mpath[,node]==1))
                                    chi <- c(tmp[which(!tmp%in%chi)], chi)
                                    isComputed[nod] <- TRUE                                    
                                }
                            }
                        }                                                
                    }
                }
                
                if(method == "TrueGraph")
                {
                    for(pa in all_parents)
                    {
                        if(interp.method == "Linear")
                        {
                            causEffect[,node] <- causEffect[,node] + approx(X[,pa], Fpa[,pa,node], xout = causEffect[,pa], method = "linear", rule = 1)$y
                        } else if(interp.method == "Spline")
                        {
                            ### USE WITH CAUTION !!! ###
                            causEffect[,node] <- causEffect[,node] + spline(X[,pa], Fpa[,pa,node], xout = causEffect[,pa])$y
                        } else
                        {
                            stop('Your specified interpolation method is not available!')
                        }
                    }
                    
                    causEffect[,node] <- causEffect[,node] + rnorm(M)*Noi[node]
                    
                } else if(method == "EstimatedGraph")
                {
                    #TT <- train_gam_causEffect(X[,all_parents],X[,node],node,TT)
                    TT[[node]] <- train_gam(X=X[,all_parents],y=X[,node])$model
                    datafr <- data.frame(causEffect[,all_parents])
                    nam <- "var2"
                    if(n_parents > 1)
                    {
                        for(ii in 2:n_parents)
                        {
                            nam <- c(nam,paste("var", ii+1, sep = ""))
                        }
                    }
                    names(datafr) <- nam
                    causEffect[,node] <- causEffect[,node] + predict.gam(object = TT[[node]], newdata = datafr , type = "link") + sample(TT[[node]]$residuals,M,replace=TRUE)
                } else if(method == "EstimatedGraphGaussianNoise")
                {
                    #TT <- train_gam_causEffect(X[,all_parents],X[,node],node,TT)
                    TT[[node]] <- train_gam(X=X[,all_parents],y=X[,node])$model
                    datafr <- data.frame(causEffect[,all_parents])
                    nam <- "var2"
                    if(n_parents > 1)
                    {
                        for(ii in 2:n_parents)
                        {
                            nam <- c(nam,paste("var", ii+1, sep = ""))
                        }
                    }
                    names(datafr) <- nam
                    causEffect[,node] <- causEffect[,node] + predict.gam(object = TT[[node]], newdata = datafr , type = "link") + rnorm(M)*sd(TT[[node]]$residuals)
                }
                
                children <- children[!children==node]
                tmp <- which(G[node,] == 1 & Mpath[,j] == TRUE)
                children <- c(tmp[which(!tmp%in%children)], children)
                isComputed[node] <- TRUE
                
            }            
        }
    }
    return(causEffect[,j])
}



find_parent_roots <- function(i,node,Mpath,G)
{
    # Find all ancestor roots of node "node"
    
    tmp <- which((Mpath[i,] == 0) & (Mpath[,node] == 1))
    roots <- which(colSums(G) == 0) 
    root_parents <- tmp[tmp %in% roots]
    
    return(root_parents)
}

