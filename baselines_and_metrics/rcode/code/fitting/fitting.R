# Copyright (c) 2010-2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
# 
# train_*** <- function(X, y, pars)
# e.g. train_linear or train_gp
#
# Performs a regression from X to y
# 
# INPUT:
#   X         nxp matrix of training inputs (n data points, p dimensions)
#   Y         vector of N training outputs (n data points)
#   pars      list containing parameters of the regression method
#
# OUTPUT:
#   result    list with the result of the regression
#      $model        list of learned model (e.g., weight vector)
#      $Yfit         fitted outputs for training inputs according to the learned model
#      $residuals    noise values (e.g., residuals in the additive noise case)



####
#Linear Regression
####
train_linear <- function(X,y,pars = list())
{
    mod <- lm(y ~ X)
    result <- list()
    result$Yfit = as.matrix(mod$fitted.values)
    result$residuals = as.matrix(mod$residuals)
    result$model = mod
    #for coefficients see list(mod$coef)
    return(result)
}



####
#GP Regression
####
gp_regression <- function(X,y, pars=list())
{
    options=gpOptions("ftc")
    options$kern$comp=list("rbf","white")
    #options$learnScales=TRUE
    model<-gpCreate(dim(X)[2],1,X,y,options)
    y2<-gpOut(model,X)
    model$Yfit<-y2
    model$residuals<-y-y2
    return(model)
}

train_gp <- function(X,y,pars = list())
{
    library(gptk)
    mod <- gp_regression(as.matrix(X),as.matrix(y))
    result <- list()
    result$Yfit = mod$Yfit
    result$residuals = mod$residuals
    result$model = mod
    return(result)
}




####
#gam Regression from mgcv package
####
library(mgcv)
train_gam <- function(X,y,pars = list(numBasisFcts = 10))
{
    if(!("numBasisFcts" %in% names(pars) ))
    { 
        pars$numBasisFcts = 10
    }
    p <- dim(as.matrix(X))
    if(p[1]/p[2] < 3*pars$numBasisFcts)
    {
        pars$numBasisFcts <- ceiling(p[1]/(3*p[2]))
        cat("changed number of basis functions to    ", pars$numBasisFcts, "    in order to have enough samples per basis function\n")
    }
    dat <- data.frame(as.matrix(y),as.matrix(X))
    coln <- rep("null",p[2]+1)
    for(i in 1:(p[2]+1))
    {
        coln[i] <- paste("var",i,sep="")
    }
    colnames(dat) <- coln
    labs<-"var1 ~ "
    if(p[2] > 1)
    {
        for(i in 2:p[2])
        {
            labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,") + ",sep="")
            #      labs<-paste(labs,"s(var",i,") + ",sep="")
            # labs<-paste(labs,"lo(var",i,") + ",sep="")
        }
    }
    labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,")",sep="")
    # labs<-paste(labs,"s(var",p[2]+1,")",sep="")
    # labs<-paste(labs,"s(var",p[2]+1,", bs = "cc")",sep="") #factor 2 faster
    # labs<-paste(labs,"s(var",p[2]+1,", bs = "cr")",sep="") # factor 2 + eps faster
    # labs<-paste(labs,"lo(var",p[2]+1,")",sep="")
    mod_gam <- FALSE
    try(mod_gam <- gam(formula=formula(labs), data=dat),silent = TRUE)
    if(typeof(mod_gam) == "logical")
    {
        cat("There was some error with gam. The smoothing parameter is set to zero.\n")
        labs<-"var1 ~ "
        if(p[2] > 1)
        {
            for(i in 2:p[2])
            {
                labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,",sp=0) + ",sep="")
            }
        }
        labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,",sp=0)",sep="")
        mod_gam <- gam(formula=formula(labs), data=dat)
    }
    result <- list()
    result$Yfit <- as.matrix(mod_gam$fitted.values)
    result$residuals <- as.matrix(mod_gam$residuals)
    result$model <- mod_gam 
    result$df <- mod_gam$df.residual     
    result$edf <- mod_gam$edf     
    result$edf1 <- mod_gam$edf1     
    
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
    return(result)
}



####
#gam boost regression from package mboost
####
library(mboost)
train_GAMboost <- function(X,y,pars = list()) #
{
    ## begin old version
    # EXPLANATION: surprisingly, it turned out that this cannot be applied to large p (private discussion with T. Hothorn in Sep 2013)
    # yy <- y
    # dat <- data.frame(cbind(yy,X))
    # gb <- gamboost(yy ~ .,data=dat, baselearner = "bbs")
    ## end old version
    
    ## begin new version
    dat <- as.data.frame(X)
    bl <- lapply(dat, bbs)
    gb <- mboost_fit(bl, y)
    ## end new version
    
    result <- list()
    result$Yfit <- gb$fitted()
    result$residuals <- gb$resid()
    result$model <- gb
    return(result)
}



####
# penGam regression 
####
train_penGAM <- function(X,y,pars = list(kCV = 3))
{
    # Meier, L., van de Geer, S. and Bühlmann, P. (2009). High-Dimensional Additive Modeling. Annals of Statistics 37, 3779-3821
    # code available from Lukas Meier [meier@stat.math.ethz.ch]
    # !! It is not so easy to find good values for cross-validation !!
    war <- library(penGAM, logical.return=TRUE)
    if(!war)
    {
        error("penGAM regression cannot be used because the package (based on Meier, L. et al. High-Dimensional Additive Modeling. 2009) cannot be loaded. Code is available from Lukas Meier [meier@stat.math.ethz.ch]")
    }
    if(is.matrix(X))
    {
        dimX <- dim(X)[1]
    } else if(is.vector(X))
    {
        dimX <- length(X)
        X <- as.matrix(X)
    }
    
    # Cross-Validation
    
    allSets <- getSetsForKFoldCV(dimX, pars$kCV)
    trainSets <- allSets$trainSets
    testSets <- allSets$testSets 
    
    lambda.pen <- 0.955^seq(1,120,by=30)
    #lambda.pen  <- 0.9^(1:20) ## 0.955^(1:110) 0.955^(1:110) ## Special!!!
    #lambda.pen  <- c(0.5,0.2,0.1) ## 0.955^(1:110) 0.955^(1:110) ## Special!!!
    #lambda.pen  <- 0.95 ## 0.955^(1:110) 0.955^(1:110) ## Special!!!
    #lambda.pen = c(1, 0.9, 0.8, 0.5, 0.3, 0.1),
    
    lambda.curv <- c(2^(4:(-4)), 0)
    #lambda.curv <- c(2^(6:(-6)), 0)
    #lambda.curv = c(6, 4, 2, 1), 
    
    SS <- matrix(0,length(lambda.pen), length(lambda.curv))
    for(i in 1:pars$kCV)
    {
        kkk <- round(sqrt(dimX))
        fit.GAM <- penGAM(as.matrix(X[trainSets[[i]],]), as.matrix(y[trainSets[[i]]]), 
                          lambda.pen,
                          lambda.curv,
                          knots = kkk, 
                          #########################
                          # is this a good idea??????
                          #########################
                          model = LinReg(), control = grpl.control(trace = 0))
        pred <- predict(fit.GAM, newdata = as.matrix(X[testSets[[i]],]))
        for(pp in 1:length(lambda.pen))
        {
            for(cc in 1:length(lambda.curv))
            {
                SS[pp,cc] <- SS[pp,cc] + sum((pred[pp,cc,] - y[testSets[[i]]])^2)
            }
        }
    }
    
    #choose best lambdas
    opt <- arrayInd(which.min(SS),dim(SS))
    lambda.pen.opt <- lambda.pen[opt[1]]
    lambda.curv.opt <- lambda.curv[opt[2]]
    cat("Optimal value for lambda.pen: ", lambda.pen.opt, "\n")
    cat("Optimal value for lambda.curv: ", lambda.curv.opt, "\n")
    
    
    # fit a gam model    
    fit.GAM <- penGAM(as.matrix(X), as.matrix(y), 
                      lambda.pen.opt,
                      lambda.curv.opt,
                      knots = round(sqrt(dim(X)[1])),
                      model = LinReg(), control = grpl.control(trace = 0))
    
    result <- list()
    result$Yfit <- predict(fit.GAM,X)
    result$residuals <- y - result$Yfit
    result$model <- fit.GAM
    return(result)
}



####
#lasso regression from package glmnet
####
war <- library(glmnet, logical.return=TRUE)
if(!war)
{
    cat("The package glmnet is not installed. This only means that lasso regression cannot be performed. However, this has no effect of the standard GAM version.\n")
}
train_lasso <- function(X,y,pars = list())
{
    cvres <- cv.glmnet(X,y)
    mod <- glmnet(X,y,lambda = cvres$lambda.1se)
    result <- list()
    result$Yfit <- predict(mod,X)
    result$residuals <- y - result$Yfit
    result$model <- mod
    return(result)
}

# =========
# 2. train_model
# =========

train_model <- function(f,X,y,pars = list())
{
    result <- f(X,y,pars)
}

