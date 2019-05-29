#This has not been yet implemented!
#maybe use d-separations and PC algorithm???


# bdag.to.mdag<-function( BDAG ) { #this would give the equivalence class MDAG
#   cat('bdag.to.mdag WARNING: NEEDS TO ADD SOME RULES FROM PEARLS BOOK!\n');
# 
#   n<-nrow(BDAG)
#   V<-vstructures(BDAG)
#   S<-skeleton(BDAG)
#   #now just create a MDAG with these properties!
#   MDAG<-array(0,dim(BDAG))
#   #first set all edges with a three in each end
#   MDAG[ S != 0 ] <- 3
# 
#   for ( j in index(1,nrow(V) ) ) { #input the v-structures!
#      #for a triple c(k,i,l)
#      k<-V[j,1]
#      i<-V[j,2]
#      l<-V[j,3]
#      #there is a directed edge from k to i
#      MDAG[i,k]<-1
#      MDAG[k,i]<-2
#      #there is a directed edge from l to i
#      MDAG[i,l]<-1
#      MDAG[l,i]<-2
#   }
# 
# 
#   oldMDAG<-array(0,dim(MDAG));
# #  A<-accessible(MDAG==1); #accessible by directed edges
# #  A2<-accessible( MDAG==3 ); #accessible by undirected edges only
#   while ( any(MDAG==3) && any(oldMDAG !=MDAG ) ) {   
#     #go through the arcs notice until no change!
#     #or no more threes
#     # print(MDAG)
#     # print(oldMDAG !=MDAG)
#     #wait(5)
#     oldMDAG<-MDAG;
# 
#     edges<-which(MDAG==3, arr.ind = TRUE)
# 
#     for ( k in 1:nrow(edges) ) {
#       b<-edges[k,1]
#       c<-edges[k,2]
#       
#       #implement rule 1 of pearls book
#       #orient b-c into b->c if there is an arrow a->b and a and c are non-adjacent      
#       for ( a in which(MDAG[b,] == 1) ) { 
#         #go through possible a:s which have an edge to b
#         if ( MDAG[a,c] == 0 ) {
#           MDAG[b,c]<-2
#           MDAG[c,b]<-1
# #          A<-accessible(MDAG==1); #update the accessible matrix
# #          A2<-accessible(MDAG==1 | MDAG==3 );
#           break;
#         }
#       }
# 
#       #implement rule 2 of pearls book
#       #orient a-b into a->b whenever there is a chain a->c->b
#       # notice that this chain could be bigger than just 3 nodes
#       a<-edges[k,1]
#       b<-edges[k,2]
# 
# #      if ( A[b,a] != 0 ) { #is this enough, chain should be only 3 nodes long!
# #        MDAG[a,b]<-2
# #        MDAG[b,a]<-1
# #        A<-accessible(MDAG==1); #update the accessible matrix
# #        A2<-accessible( MDAG==3 );
# #      }
#       for ( c in 1:n ) {
#         if (MDAG[c,a] == 1 & MDAG[b,c] == 1) {
#           MDAG[a,b]<-2
#           MDAG[b,a]<-1
# #        A<-accessible(MDAG==1); #update the accessible matrix
# #        A2<-accessible( MDAG==3 );
#           break;
#         }
#       }
# 
# #ONCE EDGE IS DIRECTED ONE SHOUD BREAK UNTIL A NEW EDGE COMES
# 
#       #the 3rd rule 
#       #orient a-b into a-> whenever there are two chains a-c->b, a-d->b
#       #such that c and d are nonadjacent
#       for ( c in 1:n ) {
#         for (d in 1:n ) {
#           if ( MDAG[c,a] == 3 &&  #there is an undirected path from a to c
#                MDAG[b,c] == 1 && #and directed b->c
#                MDAG[d,a] == 3 && #undirected a-d
#                MDAG[b,d] != 1 && #directed d->b
#                MDAG[c,d] == 0 && #c and d nonadjacent
#                MDAG[d,c] == 0 ) {
#             MDAG[a,b]<-2
#             MDAG[b,a]<-1
#             break; #break break break
#           }
#         }
#       }
#       #rule 4
#       # orient a-b into a->b when a-c->d, c->d->b, c and b nonadjacent, a nd d adjacent
#       for ( c in 1:n ) {
#         for (d in 1:n ) {
#           if ( 
# 
# 
#               ) {
#             MDAG[a,b]<-2
#             MDAG[b,a]<-1
#             break; #break break break
#           }
#         }
#       }
# 
# 
# 
#     }
# 
#     #i is the first node
#     #j is the other
#     #notice that in future rounds the i and j will be the other way around 
# 
#   }
# 
# 
#   #mdag(MDAG)
#   MDAG
# }