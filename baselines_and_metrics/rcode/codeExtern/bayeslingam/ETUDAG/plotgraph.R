#subfunctions to be used with the plotting of DAGs

plotgraph <- function( B, ngvec=rep(F,ncol(B)), filename='tempfile' ) {

  # generate the dot file
  dotfilename <- paste(sep='',filename,'.dot')
  psfilename <- paste(sep='',filename,'.ps')
  nodeshapevec <- ngvec
  if (all((B==0) | (B==1))) {
    usearclabels <- FALSE
  }
  else {
    usearclabels <- TRUE
  }

  graph2dot( B, dotfilename, nodeshapevec, usearclabels )

  # make a ps file using graphviz
  usecirco <- FALSE
  if (usecirco) useprog <- 'circo' else useprog <- 'dot'
  progcall <- paste(sep='',useprog,' -Tps ',dotfilename,' -o ',psfilename)
  system(progcall)

  # show the ps file using gv
  gvcall <- paste(sep='','gv ',psfilename)
  system(gvcall,wait=FALSE)
  wait(1.5)

}

plotgraph2start <- function( filename='tempfile' ) {
  # generate the dot file
  dotfilename <- paste(sep='',filename,'.dot')
  psfilename <- paste(sep='',filename,'.ps')

  graph2dot2start( dotfilename )

}

plotgraph2 <- function( B, ngvec=rep(F,ncol(B)), filename='tempfile', title='', comment='',nodenames=NULL ) {

  # generate the dot file
  dotfilename <- paste(sep='',filename,'.dot')
  psfilename <- paste(sep='',filename,'.ps')
  nodeshapevec <- ngvec
  if (all((B==0) | (B==1))) {
    usearclabels <- FALSE
  }
  else {
    usearclabels <- TRUE
  }

  graph2dot2( B, dotfilename, nodeshapevec, usearclabels, title, comment,nodenames=nodenames )

}

plotgraph2stop <- function( filename='tempfile' ) {
  # generate the dot file
  dotfilename <- paste(sep='',filename,'.dot')
  psfilename <- paste(sep='',filename,'.ps')

  graph2dot2stop( dotfilename )

  # make a ps file using graphviz
  usecirco <- FALSE
  if (usecirco) useprog <- 'circo' else useprog <- 'dot'
  progcall <- paste(sep='',useprog,' -Tps ',dotfilename,' -o ',psfilename)
  system(progcall)

  # show the ps file using gv
  gvcall <- paste(sep='','gv ',psfilename)
  system(gvcall,wait=FALSE)
  wait(1.5)

}


wait <- function( nseconds ) {
  i = 0
  start <- proc.time()[3]
  while (proc.time()[3] < (start+nseconds)) { i <- i+1 }

}
