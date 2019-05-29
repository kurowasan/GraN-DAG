#iperm returns the inverse order of the given order p

iperm<-function( p ) {
  q <- rep(0,length(p))
  for (i in 1:length(p)) {
    ind <- which(p==i)
    q[i] <- ind[1]
  }
  q
}