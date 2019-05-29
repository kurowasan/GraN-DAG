#changes from CDAG format to BDAG format, one or several CDAGs
bdag.to.cdag<-function( BDAGs ) {
  d<-dim(BDAGs)
  if( length(d) == 2 ) { #incase of a single BDAG handle the dimensions
    dim(BDAGs)<-c(d,1)
    d<-dim(BDAGs)
  }
  n<-d[3]
  nvars<-d[1]
  CDAGs<-array(NA,c(n,nvars+(nvars*(nvars-1)/2)))

  for (i in 1:n) {
    BDAG<-BDAGs[,,i]
    order<-causalOrder(BDAG)
    #print(order)

    BDAG<-BDAG[order,order] #now this is lower diagonal
    #print(BDAG)
    CDAGs[i,]<-c(order,BDAG[lower.tri(BDAG,diag=FALSE)])
  }
  #attr(CDAG,'class')<-'CDAG'
  CDAGs
}