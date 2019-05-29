computeScoreMatDavid <- function(X,numCores=1,parsScore=list(numBasisFcts=10),
                                 diagonalscore=0)
  # 
  #
  #
  # 2013, David Buerge, buerged@student.ethz.ch
  
{
  p<-dim(X)[2]
  
  # old (slow and not parallel)
  #
  # scorem<-matrix(NA,p,p);diag(scorem)<-0
  #
  # for(i in 1:p){
  #   for(j in 1:p){
  #     s0<--log(sum(X[,j]^2))
  #     if(i!=j)
  #     {
  #        scorem[i,j]<--log(sum((train_gam(X[,i],X[,j])$residuals)^2))-s0
  #     }
  #   }
  # }
  
  # old calculation (mean square error instead of variance)
  # score<-function(X,i,j){-log(sum((train_gam(X[,i],X[,j])$residuals)^2))+log(sum(X[,j]^2))}
  
  score<-function(X,i,j){log(var(X[,j])/var(train_gam(X=X[,i],y=X[,j],pars=parsScore)$residuals))}
  
  scorem<-mcmapply(FUN=score,MoreArgs=list(X=X),i=rep(1:p,p),j=rep(1:p,each=p),mc.cores=numCores)
  scorem<-matrix(scorem,p,p)
  
  diag(scorem)<-diagonalscore
  
  return(scorem)
  
}