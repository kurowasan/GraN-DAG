# function for plotting a hidden-variable-dag
# observed variables get round node, latent variables squared ones

plothdags <- function(hvdag) {

  n <- varnumber_hdag(hvdag)
  ngvec <- c(rep(F,n[1]),rep(T,n[2]))
  plotdags(hvdag,ngvec=ngvec,hidden=TRUE)

}