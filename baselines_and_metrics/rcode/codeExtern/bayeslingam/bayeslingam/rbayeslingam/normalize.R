# normalize:
# subtracts the mean of each column and divides by its standard deviation,
# yielding zero-mean, unit-variance variables

normalize<-function(D) {
  for ( i in 1:ncol(D)) {
    m<-mean(D[,i])
    s<-sd(D[,i])
    D[,i]<-(D[,i]-m)/s
  }
  D
}
