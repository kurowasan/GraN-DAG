#skeleton returns the skeleton of a BDAG or a MDAG

skeleton<-function( MBDAG ) {
  #DAG can be a BDAG or a MDAG
  S<-MBDAG + t(MBDAG) #make symmetric
  #attr(S,'class')<-NULL #delete the class

  S<-(S != 0)
  mode(S)<-'numeric'
  S
}