#mdag.to.bdag takes an equivalence class mdag and returns all the bdags
#in that particular equivalence class


#these result in a array of dags
#this is considerably harder task!

#so we have to input for the places of 3 ones and twos so that no new V-structures
#are generated
#one way of doing this would be to go through all dags and checking for equivalence
#this is however totally too slow for example a DAG with 20 variables
#anyway there is probably only a limited number of dags in the equivalence class!
#for 3 there are 25 dags and the largest equivalence class contains 6 DAGs
#notice that the 3 type edges can be directed so that
# 1. No new V-structures are formed
# 2. No cycles are formed

mdag.to.bdag<-function( MDAG ) { #
  #bdag( orientedges( MDAG ) )
  R<-orientedges(MDAG) == 1
  mode(R)<-'numeric'
  R
}

#orients one edge and sends the recursion one back
#returns an array of DAGs
orientedges<-function( MDAG ) {
  threes<- ( MDAG == 3 )
  if ( any( threes ) ) {
    n<-nrow(MDAG) #get the number of nodes

    nodes<-which( threes, arr.ind=TRUE ) #get the nodes that will be directed
    i<-nodes[1,1]
    j<-nodes[1,2]

    MDAG1<-MDAG2<-MDAG

    #orient the edges
    #for MDAG 1 i<-j
    MDAG1[i,j]<-1
    MDAG1[j,i]<-2

    #for MDAG 2 j<-i
    MDAG2[i,j]<-2
    MDAG2[j,i]<-1

    #originally both are valid
    MDAG1valid<-TRUE
    MDAG2valid<-TRUE

    #now make sure no new v-structures
    Vi<-vstructures(MDAG,i) #get the original v structures
    Vj<-vstructures(MDAG,j)

    Vi2<-vstructures(MDAG1,i)
    Vj2<-vstructures(MDAG2,j)

    #valid if no new v-structures!

    MDAG1valid<- ( nrow(Vi) == nrow(Vi2) && all( Vi == Vi2 ) )
    MDAG2valid<- ( nrow(Vj) == nrow(Vj2) && all( Vj == Vj2 ) )


    #check for directed cycles!
    #casting to BDAG, ignoring the threes
    if ( MDAG1valid ) MDAG1valid<- !cycles( MDAG1 == 1 )
    if ( MDAG2valid ) MDAG2valid<- !cycles( MDAG2 == 1 )


    #and no new directed cycles needed for the three node but not here


    if (MDAG1valid) { 
      MDAGs1<-orientedges( MDAG1 )
    } else {
      MDAGs1<-array(0,c(n,n,0))
    }
    if ( MDAG2valid ) {
      MDAGs2<-orientedges( MDAG2 )
    } else {
      MDAGs2<-array(0,c(n,n,0))
    }
    #combine for the end result
    R<-c( MDAGs1, MDAGs2 )

    #correct the dimensions
    dim(R)<-c(n,n,length(R)/(n^2))
    #return R
    R
  } else { #no more edges to direct
           #return just the MDAG that is a BDAG
    MDAG
  }
}