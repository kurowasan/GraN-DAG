# Copyright (c) 2010 - 2012  Jonas Peters  [jonas.peters@tuebingen.mpg.de]
# All rights reserved.  See the file COPYING for license terms. 
#
# Explanation:
#
# M: nxp matrix with n being sample size and p the number of variables 
# alpha: significance level of the independence test
# model: assumed model for regression (e.g. train_linear, train_gam or train_gp). See fitting_ts.R. 
# indtest: the independence test that should be performed (e.g. indtestHsic).
# confounder_check (only for DAG mode): if TRUE, partial causal discovery method is applied.
# check_ind_of_res: if TRUE, the method additionally checks the independence of residuals. 



fit_and_test_independence <- function(x,y,z,alpha,model,parsModel = list(),indtest, parsIndtest)
    # fits x using y and tests against z
{
    y <- as.matrix(y)
    z <- as.matrix(z)
    
    ####
    # fit x using y
    ####
    mod_fit <- train_model(model,y,x,parsModel)
    r2 <- mod_fit$resid
    
    Tquan <- indtestAll(indtest,z,r2,alpha,parsIndtest)
    return(Tquan)
}


list_subsets <- function(n,k)
{
    a<-matrix(0,n^k,k)
    for(i in 1:k)
    {
        a[,i] <- rep(1:n,each=n^(k-i),times=n^(i-1))
    }
    for(i in 1:(n^k))
    {
        if(length(unique(a[i,])) < k)
        {
            a <- a[-c(i),]
        }
    }
    a<-a
}



ICML <- function(M, alpha = 0.05, model = train_linear, parsModel = list(), indtest = indtestHsic, parsIndtest = list(method = "ExactFastTrace"), confounder_check = 0, output = FALSE)
{
    #M contains the data (each col one component)
    #confounder_check indicates subsets of which size the method tries to omit if it doesn't find any possible sink node
    stopping <- 1
    p <- dim(M)[2]
    C <- matrix(0,p,p)
    err <- matrix(0,p,1)
    S <- 1:p
    par <- matrix(0,p-1,p-1)
    parlen <- rep(0,p-1)
    variable <- rep(0,p-1)
    indtest_at_end <- rep(0,p-1)
    d <- 0
    while(length(S)>1)
    {
        #show(variable)
        d <- d+1
        # check is a vector. the k-th entry < 0 says that making the k-th variable to a sink node leads to independent residuals. 
        check <- rep(0,length(S))
        #show(S)
        for(k in 1:length(S))
        {
            i <- S[k]
            S_new <- S
            S_new <- S_new[-c(k)]
            if(output)
            {
                print(paste("fit",i,"with the help of", paste(S_new, collapse=" "),"..."))  
            }
            Fc <- fit_and_test_independence(M[,i],M[,S_new],M[,S_new],alpha,model,parsModel,indtest, parsIndtest)
            #check[k] <- Fc$statistic - Fc$crit.value
            check[k] <- -Fc$p.value
            if(output)
            {
                if(check[k]>-alpha)
                {
                #    print(paste("Independence rejected: test statistic - critical value =",check[k]))
                    print(paste("Independence rejected: p-value =",-check[k]))
                }
                else
                {
                #    print(paste("Independence not rejected: test statistic - critical value =",check[k]))
                    print(paste("Independence not rejected: p-value =",-check[k]))
                }
            }
        }
        #         if(1==0) 
        #             #        if(sum(check<0)==0) #no possible sink node found
        #         {
        #             if(confounder_check>0 && length(S)>2)
        #             {
        #                 show("Since no possible sink node was found, the algorithm tries to omit dimensions...")
        #                 for(sizesubset in 1:confounder_check)
        #                 {
        #                     print(paste("tries to omit", sizesubset, "dimension(s) of the time series..."))
        #                     a <- list_subsets(length(S),sizesubset)
        #                     pp <- dim(a)		
        #                     check2 <- matrix(0,pp[1],length(S)-pp[2])
        #                     show("Does omitting variables help?")
        #                     for(k in 1:pp[1])
        #                     {
        #                         #S[a[k,]] werden entfernt
        #                         S_new <- S
        #                         S_new <- S_new[-a[k,]]
        #                         for(kk in 1:length(S_new))
        #                         {
        #                             #Is i possible sink?
        #                             i <- S_new[kk]	                
        #                             S2_new <- S_new
        #                             S2_new <- S2_new[-c(kk)]
        #                             Fc <- fit_and_test_independence(M[,i],M[,S2_new],M[,S2_new],alpha,model,parsModel,indtest,parsIndtest)
        #                             check2[k,kk] <- Fc[1]-Fc[2]
        #                             print(paste("omitting: ", paste(S[a[k,]],collapse=" "), "Sink ", i, " leads to ", Fc[1]-Fc[2], " (<0 independence)."))
        #                         }
        #                     }
        #                     if(sum(sum(check2<0))>0)
        #                     {
        #                         #found something!
        #                         stopping <- 0
        #                         pp <- dim(check2)
        #                         k1 <- which.min(check2)
        #                         k2 <- (k1-1)%/%pp[1]+1
        #                         k1 <- k1%%pp[1]
        #                         if(k1==0)
        #                         {
        #                             k1 <- pp[1]
        #                         }
        #                         S_new <- S
        #                         S_new <- S_new[-a[k1,]]
        #                         variable[(d):(d+sizesubset-1)] <- S[a[k1,]]
        #                         err[(d):(d+sizesubset-1)] <- rep(1,sizesubset)
        #                         for(iii in 1:sizesubset)
        #                         {
        #                             for(jjj in 1:length(S))
        #                             {
        #                                 C[S[a[k1,iii]],S[jjj]] <- -1
        #                                 C[S[jjj],S[a[k1,iii]]] <- -1
        #                             }
        #                         }
        #                         
        #                         variable[d+sizesubset] <- S_new[k2]
        #                         S <- S_new[-k2]
        #                         parlen[d+sizesubset] <- length(S)
        #                         par[d+sizesubset,1:length(S)] <- S
        #                         
        #                         d<-d+sizesubset
        #                         break
        #                     }	    
        #                 } # end for
        #                 
        #                 if(stopping==1)
        #                 {
        #                     show("Not even omitting variables helped. Stop the search.")
        #                     err[d:(p-1)] <- rep(1,(p-d))
        #                     for(i in 2:length(S))
        #                     {
        #                         for(j in 1:(i-1))
        #                         {
        #                             C[S[i],S[j]] <- -1
        #                             C[S[j],S[i]] <- -1
        #                         }
        #                     }
        #                     break
        #                 }
        #                 check2 <- rep(0,p-1)
        #                 stopping <- 1
        #             }
        #             else # no possible sink node and confounder_check disabled
        #             {
        #                 show("No possible sink node found. Stop the search.")
        #                 err[d:(p-1)] <- rep(1,(p-d))
        #                 for(i in 2:length(S))
        #                 {
        #                     for(j in 1:(i-1))
        #                     {
        #                         C[S[i],S[j]]<-88
        #                         C[S[j],S[i]]<-88
        #                     }
        #                 }
        #                 break
        #             }
        #         }
        #         else #possible sink node found
        #         {
        bb <- which.min(check)
        variable[d] <- S[bb]
        S <- S[-c(bb)]
        parlen[d] <- length(S)
        par[d,1:length(S)] <- S
        if(output)
        {
            print(paste("Possible sink node found:",variable[d]))
            print(paste("causal order (beginning at sink):",paste(variable,collapse=" ")))
        }	
        #         }
        rm(check)
    }
    # show(variable)
    if(d<p)
    {
        variable[p]<-S[1]
    }
    if(output)
    {   
        print(paste("causal order (beginning at sink):",paste(variable,collapse=" ")))
        print(paste("removing unnecessary edges..."))
    }
    rm(S)
    #todo: here, we take the first possible parent away (not the best one). in theory it probably doesn't make a difference. experiments in paper done correctly. but maybe finding all dags is better.
    for(d in 1:(p-1))
    {
        if(err[d] != 1)
        {
            S<-par[d,1:parlen[d]]
            for(i in 1:length(S))
            {
                S_new<-S
                S_new<-S_new[-c(1)]
                if(length(S)==1)
                {
                    tsx <- M[,variable[d]]
                    # todo: test against par[d,1:parlen[d]] or S????
                    Fc <- indtestAll(indtest,M[,par[d,1:parlen[d]]],tsx,alpha,parsIndtest)
                }
                else
                {
                    # todo: test against par[d,1:parlen[d]] or S????
                    Fc <- fit_and_test_independence(M[,variable[d]],M[,S_new],M[,par[d,1:parlen[d]]],alpha,model,parsModel,indtest, parsIndtest)
                }
                if(Fc$statistic < Fc$crit.value)
                {
                    S <- S_new
                }
                else
                {
                    if(length(S)>1)
                    {
                        tmp <- S[1]
                        S[1:(length(S)-1)] <- S[2:length(S)]
                        S[length(S)] <- tmp
                    }
                }
            }
            
            #todo: always hsic?
            #todo: change code...
            if(1==0)
            {                 
                print(paste("...and performing final independence test using HSIC..."))
                if(length(S)>0)
                    #if(indtest != indtestts_hsic & length(S)>0)
                {
                    print(paste("fitting ", variable[d], " with the help of ", paste(S,collapse=" "), " and testing independence against ", paste(par[d,1:parlen[d]],collapse=" ")))
                    Fc <- fit_and_test_independence(M[1:p[1],variable[d]],M[1:p[1],S],M[1:p[1],par[d,1:parlen[d]]],alpha/(p-1),max_lag,model,indtest,instant)
                }
                else
                {
                    print(paste("fitting ", variable[d], " with the help of NOTHING and testing independence against ", paste(par[d,1:parlen[d]],collapse=" ")))
                    tsx <- M[1:p[1],variable[d]]
                    pars <- list(maxOrder = max_lag)
                    modx <- traints_model(model,tsx,list(),pars)
                    resx <- modx$residuals[(modx$model$order+1):length(tsx)]
                    Fc <- indtest_model(indtest,M[1:p[1],par[d,1:parlen[d]]],resx,alpha,max_lag,FALSE)
                }
                show(sprintf("Test statistic: %.3f and critical value: %.3f and p-value %.2e", Fc[1],Fc[2],Fc[3]))        
                indtest_at_end[d] <- sign(Fc[1]-Fc[2])
            }
            else
            {
                #-1: ind., +1 dep.
                indtest_at_end[d] <- -1
            }
            parlen[d] <- length(S)
            C[S,variable[d]] <- rep(1,length(S))
        }
        else
        {
            #-1: ind., +1 dep.
            indtest_at_end[d] <- -1
        }
    }
    if(max(indtest_at_end)<0)
    {
        if(output)
        {
            print(paste("all correct..."))
        }
        return(C)
    }
    else
    {
        print(paste("final ind. test failed. No solution."))
        return(NULL)
    }
}


