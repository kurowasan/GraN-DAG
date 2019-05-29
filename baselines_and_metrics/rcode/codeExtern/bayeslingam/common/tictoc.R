#matlab tictoc functionality

tictic<-0

tic<-function() {
  tictic<<-proc.time()[3]
}

toc<-function() {
  proc.time()[3]-tictic
}