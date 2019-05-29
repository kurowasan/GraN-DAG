#converts a tdag "[X1][X2|X1]" to a bdag

tdag.to.bdag<-function(TDAGs) {
  if ( is.null(dim(TDAGs)) ) {
    dim(TDAGs)<-c(1,1)
  }
  nvars<-NULL
  B<-NULL
  varnames<-NULL
  #print(TDAGs)
  #print(B)
  #print(dim(B))

  for ( kk in 1:nrow(TDAGs) ) {
    TDAG<-TDAGs[kk,] 
   #the names of the variables in order!
    s<-substr(TDAG,start=2,stop=nchar(TDAG)-1) 
    #first rip off the beginning [ and ending ]
    s<-strsplit(s,"\\]\\[")[[1]]#the split with ][
    #print(s)
    if ( is.null(nvars) ) {
      nvars<-length(s)
      B<-array(0,c(nvars,nvars,nrow(TDAGs)))
      varnames<-paste("X",c(1:nvars),sep="")
    }
    for (i in 1:length(s)) {
      si<-s[i]
      k<-strsplit(si,"\\|")[[1]]
      node<-which(varnames==k[1])
      if ( length(k) > 1 ) {
        parentsstring<-k[2]
      } else {
        parentsstring<-""
      }
      parents<-c()
      while( nchar(parentsstring) != 0 ) {
        for (j in 1:length(varnames)) {
          if ( substr(parentsstring,start=1,
                      stop=nchar(varnames[j])) == varnames[j] ) {
            parents<-c(parents,j)
            parentsstring<-substr(parentsstring,start=nchar(varnames[j])+2,
                                  stop=nchar(parentsstring))
          }
        }
      }
      #print(B)
      #cat('node,parents,k:');print(c(node,parents,kk))
      B[node,parents,kk]<-1
    }
  }
  B
}