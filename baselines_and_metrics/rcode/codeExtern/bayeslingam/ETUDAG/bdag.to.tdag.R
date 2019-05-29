#takes BDAGs and transforms them to text based tdags

bdag.to.tdag<-function(BDAGs) {
  d<-dim(BDAGs)
  if( length(d) == 2 ) { #incase of a single BDAG handle the dimensions
    dim(BDAGs)<-c(d,1)
    d<-dim(BDAGs)
  }
  n<-d[3]
  nvars<-d[1]
  TDAGs<-rep("",n)
  dim(TDAGs)<-c(n,1)

  for ( i in 1:n) {
    for ( node in 1:nvars) {
      pa<-which(BDAGs[node,,i]==1)
      if ( !is.null(pa) && length(pa) > 0 ) {
        pastring<-paste(sprintf(':X%i',pa),collapse='')
        substr(pastring,1,1)<-'|'
      } else {
        pastring<-''
      }
      #cat('pastring:');print(pastring)
      #cat('c:');print(c(TDAGs[i,],'[X',node,pastring,']',sep=''))
      #cat('paste:');print(paste(TDAGs[i,],'[X',node,pastring,']',sep=''))

      TDAGs[i,]<-paste(TDAGs[i,],'[X',node,pastring,']',sep='')
    }
  }
  TDAGs
}