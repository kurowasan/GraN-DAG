#changes from CDAG format to BDAG format, one or several CDAGs

cdag.to.bdag<-function( CDAGs ) {

  if ( is.vector(CDAGs) ) {
    dim(CDAGs)<-c(1,length(CDAGs))
  }

  nvars<-max(CDAGs)
  n<-nrow(CDAGs)

  BDAGs<-vector('numeric',nvars*nvars*n) #might want this as 'logical'
  dim(BDAGs)<-c(nvars,nvars,n)


  for ( i in 1:n ) {
    #B<-vector('logical',nvars*nvars) #logical for large storing!
    B<-array(0,c(nvars,nvars))

    B[lower.tri(B,diag=FALSE)] <- CDAGs[i,nvars+(1:(nvars*(nvars-1)/2))]

    order<-CDAGs[i,1:nvars]
    invorder<-iperm(order)
    B<-B[invorder,invorder]
    #attr(B,'class')<-'BDAG'
    #mode(B)<-'logical' #for storing!
    BDAGs[,,i]<-B
  }
  if ( n == 1 ) {
    dim(BDAGs)<-dim(BDAGs)[1:2]
  }
  BDAGs
}