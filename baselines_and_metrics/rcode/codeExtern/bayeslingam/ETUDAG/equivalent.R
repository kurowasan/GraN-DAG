#whether two objects belong to the same equivalence class
# the object might be mdags or bdags
# test for skeleton and v-structures, so strictly speaking
# no extra constraints to an equivalence class are allowed for the mdag, mdag case

equivalent<-function( MBDAG1, MBDAG2 ) {
  #the dags can be MDAGS or BDAGS!
  R<-FALSE
  if ( nrow(MBDAG1) == nrow(MBDAG2) ) { #must have same dimensions!
    if ( all(skeleton(MBDAG1) == skeleton(MBDAG2) ) ) { #same skeletons
      V1<-vstructures(MBDAG1)
      V2<-vstructures(MBDAG2)
      if ( nrow(V1) == nrow(V2) ) { #same number of v structures
        if ( all ( V1 == V2 ) ) { #same v structures!
          R<-TRUE
        }
      }
    } 
  }
  R
}

