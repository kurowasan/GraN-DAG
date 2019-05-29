#test whether a give bdag candidate has cycles

cycles<-function( BDAG ) { #with MDAG this is not needed
  A<-BDAG #
  #cat('Testing for cycles:\n')
  #print(BDAG)
  Aprev<-array(0,dim(BDAG))
  C<-any(diag(A) != 0 )
  while( !C & any( A != Aprev ) ) {
    Aprev<-A
    for ( i in 1:nrow(A) ) {
      Ai<-Aprev[ ,c(i,which(Aprev[,i]==1))]
      if ( is.array(Ai) ) { #if many cols, or them up
        Ai<-apply(Ai,1,sum) != 0
      }
      A[,i]<-Ai
    }
    #print(A)
    C<-any(diag(A) != 0 )
  }
  #print(C)
  C
}

