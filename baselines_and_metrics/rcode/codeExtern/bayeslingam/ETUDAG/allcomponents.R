# allcomponents:
# returns all components (families) for a directed graph in a list
# where the node is the node integer number, and the parents are
# denoted by '1':s whereas non-parents are denoted by '0':s.
#
# INPUT:
# N - number of nodes
#

allcomponents <- function(N) {

  # create an array which holds all 2^N binary vectors
  edgeconfs<-array(NA,c(2^N,N)) 
  for ( i in 1:(2^N)) {
    edgeconfs[i,]<-dec.to.bin(i-1,N)
  }

  # this is the number of parent configurations for any single node
  k <- 2^(N-1) 

  # initializes the nodes and edges arrays
  nodes<-rep(0,N*k) 
  edges<-array(NA,c(N*k,N)) 

  # go through each node
  for (i in 1:N) { 

    # set these k elements to i 
    index <- ((i-1)*k+1):(i*k) 
    nodes[index] <- i

    # only those configurations are valid where no X->X
    index2<-which(edgeconfs[,i] == 0 ) 
    edges[index,]<-edgeconfs[index2,] 
  }

  # return a list with the nodes and edges
  list(node=nodes,edges=edges)
}
