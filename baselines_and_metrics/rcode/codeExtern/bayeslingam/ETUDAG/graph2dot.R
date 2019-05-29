graph2dot <- function( adj, filename, nodeshapevec, usearclabels ) {

  # number of nodes in graph
  nnodes <- dim(adj)[1]

  # node names currently just 1:nnodes
  nodelabels <- list()
  for (i in 1:nnodes) {
    nodelabels[[i]] <- sprintf("%d",i)
  }
    
  # node shapes: boolean vector indicating FALSE=round or TRUE=square 
  nodeshapes <- list()
  for (i in 1:nnodes) {
    if (nodeshapevec[i]) {
      nodeshapes[[i]] <- 'box'
    }
    else {
      nodeshapes[[i]] <- 'ellipse'
    }
  }
  
  # usearclabels is boolean, if TRUE we put numerical values on arcs
  if (usearclabels) {
    arclabels <- list()
    for (i in 1:(nnodes*nnodes)) {
      arclabels[[i]] <- sprintf("%.3f",adj[i])
    }
  }

  # construct a format string for nodes
  nodeformat <- '  %d [label = "%s", shape = %s, style = "%s"];'

  # construct a format string for edges
  if (usearclabels) {
    attributes <- 'style = %s, color = %s, label = "%s", dir = %s'
  }
  else {
    attributes <- 'style = %s, color = %s, dir = %s'
  }
  edgeformat = paste(sep='','  %d -> %d [', attributes, '];')

  # if filename is empty then write to a tempfile
  if (filename == "") {
    filename <- "tempfile.dot"
  }
  
  # write the beginning of the file
  write('digraph G {', file=filename, append=FALSE, sep="")
  write('  center = 1;', file=filename, append=TRUE, sep="")
  write('  size = \"10, 10\";', file=filename, append=TRUE, sep="")

  # write the nodes
  for (node in 1:nnodes) {
    thisnode <- sprintf(nodeformat, node, nodelabels[[node]],
                        nodeshapes[[node]], 'solid')
    write(thisnode, file=filename, append=TRUE, sep="")
  }

  
  # write the edges
  for (node1 in 1:nnodes) {
    arcs <- which(adj[node1,]!=0)
    for (node2 in arcs) {
      style <- 'solid'
      color <- 'black'
      dowrite <- TRUE
      if (usearclabels) {
        thisarc <- sprintf(edgeformat, node2, node1, style, color,
                           arclabels[[(node2-1)*nnodes + node1]], 'forward')
      }
      else {
        if (adj[node2,node1]!=0) {
          if (node2<node1) {
            color <- 'red'
            thisarc <- sprintf(edgeformat, node2, node1, style, color, 'none' )
          }
          else
            dowrite <- FALSE
        }
        else
          thisarc <- sprintf(edgeformat, node2, node1, style, color, 'forward')
      }
      if (dowrite) {
        write(thisarc, file=filename, append=TRUE, sep="")
      }
    }
  }
  write('}', file=filename, append=TRUE, sep="")
}

graph2dot2start <- function( filename ) {
  #print('starting')
  #print(filename)

  # write the beginning of the file
  write('digraph G {', file=filename, append=FALSE, sep="")
  write('  center = 1;', file=filename, append=TRUE, sep="")
  write('  size = \"10, 10\";', file=filename, append=TRUE, sep="")


}

graph2dot2 <- function( adj, filename, nodeshapevec, usearclabels, title, comment,nodenames=NULL ) {

  # number of nodes in graph
  nnodes <- dim(adj)[1]

  # node names currently just 1:nnodes
  nodelabels <- list()
  for (i in 1:nnodes) {
    nodelabels[[i]] <- sprintf("%s_%d",title,i)
  }
  # node shapes: boolean vector indicating FALSE=round or TRUE=square 
  nodeshapes <- list()
  for (i in 1:nnodes) {
    if (nodeshapevec[i]) {
      nodeshapes[[i]] <- 'box'
    }
    else {
      nodeshapes[[i]] <- 'ellipse'
    }
  }
  
  # usearclabels is boolean, if TRUE we put numerical values on arcs
  if (usearclabels) {
    arclabels <- list()
    for (i in 1:(nnodes*nnodes)) {
      arclabels[[i]] <- sprintf("%.3f",adj[i])
    }
  }

  # construct a format string for nodes
  nodeformat <- '  %s [label = "%s", shape = %s, style = "%s"];'

  # construct a format string for edges
  if (usearclabels) {
    attributes <- 'style = %s, color = %s, label = "%s", dir = %s'
  }
  else {
    attributes <- 'style = %s, color = %s, dir = %s'
  }
  edgeformat = paste(sep='','  %s -> %s [', attributes, '];')

  # if filename is empty then write to a tempfile
  if (filename == "") {
    filename <- "tempfile.dot"
  }
  #print('writing')
  #print(filename)
  write(sprintf('  subgraph cluster%s {',title), file=filename,append=TRUE, sep="")

  write(sprintf('  label = \"%s\"',comment), file=filename,append=TRUE, sep="")

  # write the nodes
  for (node in 1:nnodes) {
    if ( is.null(nodenames) ) {
      thisnode <- sprintf(nodeformat, nodelabels[[node]], node,
                        nodeshapes[[node]], 'solid')
      write(thisnode, file=filename, append=TRUE, sep="")
    } else {
      thisnode <- sprintf(nodeformat, nodelabels[[node]], nodenames[node],
                        nodeshapes[[node]], 'solid')
      write(thisnode, file=filename, append=TRUE, sep="")
    }

  }

  
  # write the edges
  for (node1 in 1:nnodes) {
    arcs <- which(adj[node1,]!=0)
    for (node2 in arcs) {
      style <- 'solid'
      color <- 'black'
      dowrite <- TRUE
      if (usearclabels) {
        thisarc <- sprintf(edgeformat, nodelabels[[node2]], nodelabels[[node1]], style, color,
                           arclabels[[(node2-1)*nnodes + node1]], 'forward')
      }
      else {
        if (adj[node2,node1]!=0) {
          if (node2<node1) {
            color <- 'red'
            thisarc <- sprintf(edgeformat, nodelabels[[node2]], nodelabels[[node1]], style, color, 'none' )
          }
          else
            dowrite <- FALSE
        }
        else
          thisarc <- sprintf(edgeformat, nodelabels[[node2]], nodelabels[[node1]], style, color, 'forward')
      }
      if (dowrite) {
        write(thisarc, file=filename, append=TRUE, sep="")
      }
    }
  }
  write('}', file=filename, append=TRUE, sep="")
}

graph2dot2stop <- function( filename ) {
  #print('stopping')
  #print(filename)
  write('}', file=filename, append=TRUE, sep="")
}

