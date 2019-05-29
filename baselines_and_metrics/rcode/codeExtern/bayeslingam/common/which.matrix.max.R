which.matrix.max<-function(Z) {
#this returns the index of the largest element for a matrix
  rowmax<-apply(Z,1,which.max)
  i<-which.max(Z[cbind(1:length(rowmax),rowmax)])
  j<-rowmax[i]

  c(i,j)
}
