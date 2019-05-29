#returns one of the causal orders of BDAG

causalOrder<-function(BDAG) {
  B<-BDAG

  #this finds the causal order
  order<-1:nrow(B)
  labels<-1:nrow(B)
  for ( i in 1:(nrow(BDAG)-1) ) {
    j<-which(rowSums(B) == 0)[1]

    pre<-order[index(1,i-1)]
    rest<-order[index(i,length(order))]

    #take the j:th out of the rest
    out<-rest[j]
    rest<-rest[-j]

    order<-c(pre,out,rest)
#    helper<-order[i]
#    order[i]<-order[j+i-1] #put the i:th variable to the place
#    order[j+i-1]<-helper

    #B<-B[-j,-j] #delete that row from B and reorder! 
    #order[(i+1):nrow(B)]

    B<-B[-j,-j]
    #print(order)
  }

  order
}