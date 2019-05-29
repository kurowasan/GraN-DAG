#accessible returns a matrix whose i:th column contains 
#the number of steps needed to access the node of the row of a BDAG

accessible<-function( BDAG ) {
  A<-BDAG #
  #mode(A)<-'integer'
  #class(A)<-NULL

  A[A==0]<-Inf
  #cat('Testing for cycles:\n')
  #print(BDAG)
  #print(A)
  Aprev<-array(0,dim(BDAG))
  while( any( A != Aprev ) ) { #while A changes
    Aprev<-A
    for ( i in 1:nrow(A) ) {
      ind<-which(is.finite(Aprev[,i])) #which nodes are accessible
      #print(ind)
      Ai<-as.matrix(Aprev[ ,ind])
      #print(Ai)

      for ( j in index(1,length(ind)) ) {
        Ai[,j]<-Ai[,j]+Aprev[ind[j],i] #add how many steps there was to get to the node
      }

      Ai<-cbind(Aprev[,i],Ai)

     #then choose the minimum, but non zero of the columns
      Ai<-apply(Ai,1,min)
      A[,i]<-Ai
    }
    #print(A)
  }
  A[is.infinite(A)]<-0
  A
}