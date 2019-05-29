#Returns all the neighbouring DAGs (1 turn, addition, deletion of an edge),
#to a BDAG, and the BDAG itself!
#Should the original DAG also be returned??
#An upper bound for the number of neighbouring dags
#> n*n-n
#  [1]    0    2    6   12   20   30   42   56   72   90  110  132  156  182  210
# [16]  240  272  306  342  380  420  462  506  552  600  650  702  756  812  870
# [31]  930  992 1056 1122 1190 1260 1332 1406 1482 1560 1640 1722 1806 1892 1980
# [46] 2070 2162 2256 2352 2450 2550 2652 2756 2862 2970 3080 3192 3306 3422 3540
# [61] 3660 3782 3906 4032 4160 4290 4422 4556 4692 4830 4970 5112 5256 5402 5550
# [76] 5700 5852 6006 6162 6320 6480 6642 6806 6972 7140 7310 7482 7656 7832 8010
# [91] 8190 8372 8556 8742 8930 9120 9312 9506 9702 9900



neighbourDAGs<-function( BDAG ) {
  
  #three ways to find neighbours
  
  #1. delete an edge = delete ones in BDAG
  
  #2. add an edge = add an 1 where there is a 0,
  # Check: no self cycles so b(i,j) cannot be 1 if b(j,i) is 1
  # Check: no cycles result, this just must be checked
  
  #3. turn  turn an edge i->j to i<-j
  #   Check: no cycles in the result
  
  #4. Do nothing, A dag is itselfs own neighbour

  #So it is better to go through all elements i,j
  #if i=j skip
  # if b(i,j)=1 delete it
  # if b(i,j)=0 and b(j,i) = 1 turn the edge
  # if b(i,j)=0 and (b(j,i) = 0 then add edge i->j
  
  #Create an 3-dimensional matrix of the b-matrices!
  nvars<-ncol(BDAG)
  
  #each element corresponds to a possible neighbour except for the diagonal ones,
  # + and original 1
  n<-(nvars*nvars-nvars+1)
  BDAGs<-vector('numeric',n*nvars*nvars) #might want this as 'logical'
  dim(BDAGs)<-c(nvars,nvars,n)
  
  BDAGs[,,1]<-BDAG #input the original!
  Bindex=2;
  Bcycles<-rep('FALSE',n)
  
  for (i in 1:nvars) {
    for (j in 1:nvars) {
      if ( i != j ) { #dont bother with the diagonal elements
        BDAGs[,,Bindex]<-BDAG
        if ( BDAG[i,j] == 1 ) { #delete the edge
          BDAGs[i,j,Bindex]<-0;
        } else if ( BDAG[j,i] == 1 ) { #turning the edge
          BDAGs[i,j,Bindex]<-1;
          BDAGs[j,i,Bindex]<-0;
          Bcycles[Bindex]<-cycles(BDAGs[,,Bindex]) #check for cycles!
        } else { #now an edge is just added 
          BDAGs[i,j,Bindex]<-1
          Bcycles[Bindex]<-cycles(BDAGs[,,Bindex]) #check for cycles!
        }
        Bindex<-Bindex+1 #go to the next one
      }
    }
  }
  #return only the ones where no cycles!!!
  BDAGs[,,Bcycles==FALSE]
}




#Should use cycles-function or some sort of test function???
#use cycles from the start at least


#how about just adding the stuff and later checking for cycles?
