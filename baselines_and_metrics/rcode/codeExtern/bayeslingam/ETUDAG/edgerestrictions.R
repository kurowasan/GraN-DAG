edgerestrictions <- function( allDAGs, noedge ) {

  # noedge is a matrix indicating between which pairs of nodes there are no edges
  # allowed, f. ex: noedge = 1 2
  #                          2 4
  #                          3 4
  #          means no edge between 1 and 2, between 2 and 4 and between 3 and 4.

  for (i in 1:nrow(noedge)) {
    bool <- allDAGs[noedge[i,1],noedge[i,2],]==0 &
            allDAGs[noedge[i,2],noedge[i,1],]==0
    allDAGs <- allDAGs[,,bool]
    if (length(dim(allDAGs)) == 2) dim(allDAGs) <- c(dim(allDAGs),1)
  }

  allDAGs

}

hiddenedgerestrictions <- function( allhvDAGs, noedge ) {

  n <- varnumber_hdag(allhvDAGs)

  for (i in 1:nrow(noedge)) {
    for (j in (n[1]+1):(n[1]+n[2])) {
      bool <- allhvDAGs[noedge[i,1],j,]==1 & allhvDAGs[noedge[i,2],j,]==1
      allhvDAGs <- allhvDAGs[,,!bool]
      if (length(dim(allhvDAGs)) == 2) dim(allhvDAGs) <- c(dim(allhvDAGs),1)
    }
  }

  allhvDAGs

}