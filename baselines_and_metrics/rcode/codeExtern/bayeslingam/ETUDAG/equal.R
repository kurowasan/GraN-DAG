#defining equal operators!
#check also dimensions!


equal<-function( BDAG1, BDAG2 ) {
#   R<-FALSE
#   if ( all(dim(BDAG1) == dim(BDAG2)) ) {
#     if ( all(BDAG1 == BDAG2) ) {
#       R<-TRUE
#     }
#   }
#   R
  all(dim(BDAG1) == dim(BDAG2)) && all(BDAG1 == BDAG2)

}

equal.MDAG<-function( MDAG1, MDAG2 ) {
#   R<-FALSE
#   if ( all(dim(MDAG1) == dim(MDAG2)) ) {
#     if ( all(MDAG1 == MDAG2) ) {
#       R<-TRUE
#     }
#   }
#   R
  all(dim(MDAG1) == dim(MDAG2)) && all(MDAG1 == MDAG2)

}

#'==.BDAG'<-function( BDAG1, BDAG2 ) { all( as.vector(BDAG1) == as.vector(BDAG2) ) }
#'!=.BDAG'<-function( BDAG1, BDAG2 ) { any( as.vector(BDAG1) != as.vector(BDAG2) ) }


#these must be first transformed into bdags actually
#

#'==.CDAG'<-function( CDAG1, CDAG2 ) { stop('Comparison of CDAG objects not valid!') }
#'!=.CDAG'<-function( CDAG1, CDAG2 ) { stop('Comparison of CDAG objects not valid!') }

