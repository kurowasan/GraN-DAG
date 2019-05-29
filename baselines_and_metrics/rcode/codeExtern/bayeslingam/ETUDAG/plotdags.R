# plotdags plots a list of dags, subtitles is a vector of subtitles for each
# DAG
# try plotdags(alldags(3),subtitles=sprintf('%i',1:25))
# ( try plotdags(cdag.to.bdag(alldags(3)),subtitles=sprintf('%i',1:25)) )

plotdags<-function(BDAGs, subtitles=NULL, ngvec=rep(F,ncol(B)), hidden=FALSE) {
  d<-dim(BDAGs)
  nvars<-d[1]
  if (length(d) == 2 ) {
    dim(BDAGs)<-c(d,1)
    d<-dim(BDAGs)
  }

  if ( is.null(subtitles)) {
    subtitles<-rep('',d[3])
  }

  # fcts "plotgraph2start", "plotgraph2", "plotgraph2stop" in file plotgraph.R
  plotgraph2start()

  for (i in (1:d[3]) ) {
  
    B <- BDAGs[,,i]
  
     # BDAGs is matrix of HDAGs (ie there are hidden variables)
     if (hidden) {
       # cut out hidden variables which aren't connected to the observed var.
       n <- varnumber_hdag(BDAGs[,,i])
       temp <- NULL
       for (j in (n[1]+1):d[2]) {
         if (all(BDAGs[,j,i]==0)) temp <- c(temp,j)
       }
       B <- BDAGs[setdiff(1:nvars,temp),setdiff(1:nvars,temp),i]
     }

     plotgraph2( B, ngvec, title=sprintf('g%g',i), comment=subtitles[i],
              nodenames=NULL )
     #CDAG<-bdag.to.cdag(BDAGs[,,i])
     #plotmydag2( CDAG, nvars, sprintf('g%g',i), 
     #           subtitles[i],nodenames=NULL,ngvec)
  }

  plotgraph2stop()
}

#-----------------------------------------------------------------------------
# plotmydag2:
# plots one DAG using graphviz
#-----------------------------------------------------------------------------

plotmydag2 <- function( thedag, nvars, title, comment,nodenames=NULL, ngvec=rep(F,max(thedag)) ) {

  ip <- iperm(thedag[1:nvars])
  B <- matrix(0,nvars,nvars)
  B[lower.tri(B,diag=FALSE)] <- thedag[nvars+(1:(nvars*(nvars-1)/2))]
  B <- B[ip,ip]
  plotgraph2( B, ngvec, title=title, comment=comment,
              nodenames=nodenames )

}
