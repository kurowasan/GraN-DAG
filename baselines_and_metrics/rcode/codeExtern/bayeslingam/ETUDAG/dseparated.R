# dseparated returns TRUE if x is d-separated from y given Z
# (Z possibly empty, indicated with Z=NULL)

dseparated <- function( MBHDAG, x, y, Z=NULL) {

  # checking if nodes x and y are d-connected given Z
  # at the moment kind of depth-first search, I think...

  # is a good way of doing this a depth first search or breadth first search?
  # where you can go from one node to another only if you came with an end,
  # or leave with an end or the node is in the given set!

  B <- MBHDAG

  x.par <- which(B[x,]!=0) # parents of x
  x.chi <- which(B[,x]!=0) # children of x

  dcon <- FALSE

  # check if there is a path to y via the parents of x
  for (j in x.par) {
    dcon <- dconnected(B, j, y, Z, FALSE, x)
    if (dcon) break # if I got one TRUE then stop for-loop over j
  }
  if (dcon) return(FALSE) # if I got one TRUE then x and y are d-connected
                          # given Z that means not d-separated given Z

  # check if there is a path to y via the children of x
  for (j in x.chi) {
    dcon <- dconnected(B, j, y, Z, TRUE, x)
    if (dcon) break # if I got one TRUE then stop for-loop over j
  }
  if (dcon) return(FALSE) # got a TRUE: x,y d-connected given Z

  # if non of the above found a dconnecting path: get d-separation:
  if (!dcon) return(TRUE)

}



