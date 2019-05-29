#vstructures returns the v-structures of a BDAG or MDAG
#parameter lines can be modified so that v-structures are only
#searched from given lines

vstructures<-function( MBDAG, lines=1:nrow(MBDAG) ) {
  #DAG can be a BDAG or a MDAG
  #lines determines on which line the v-structure is searched in
  #assumes 

  #returns the v-structures in a array
  #where (k,i,l) means k->i<-l and k < l

  Vs<-array(0,c(0,3))

  #lines 1:nrow(MBDAG)
  for ( i in  lines ) { #for each variable
#    cat('node:',i)
    #these are the variables that there are edges from

    #for the MDAG this must be 1, 2 or 3 is not ok
    pa<- which(MBDAG[i,] == 1) 
    #cat('parents:',pa)
    for ( k in index(1,length(pa)) ) { #k is the first node
      for (l in index(k+1,length(pa)) )  { #l is the second node

         #now if there is not a edge from k to l or l to k
         # this is a v-structure!
       # cat('considering:',pa[k],i,pa[l])
        if ( MBDAG[pa[l],pa[k]] == 0 & MBDAG[pa[k],pa[l]] == 0 ) { #it is a v-structure!
          Vs<-rbind(Vs,c(pa[k],i,pa[l]))
        }#if v
      }#for l
    }#for k
  }#for i
  Vs

}



