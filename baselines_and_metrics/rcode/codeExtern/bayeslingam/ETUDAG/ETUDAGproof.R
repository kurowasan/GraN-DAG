#function to test how the package functions work, and possibly prove that
#they work correctly

ETUDAGproof<-function() {
#   cat('Creating all 5 node dags: A<-alldags(5)\n');
#   A<-alldags(5)
# 
#   cat('Changing to BDAG format: B<-cdag.to.bdag(A)\n');
#   B<-cdag.to.bdag(A)
# 
#   cat('Changing back to CDAG format: C<-bdag.to.cdag(B)\n')
#   C<-bdag.to.cdag(B)
# 
#   cat('All same: all( C == A)\n')
#   print(all(C == A))
#   A<-alldags(3)
#   for (i in 1:25 ) {
#     print(A[i,])
#     b<-cdag.to.bdag(A[i,])
#     #print(b)
#     m<-bdag.to.mdag(b)
#     #print(m)
#     bs<-mdag.to.bdag(m)
#     #print(bs)
#     cat(' is equivalent to:\n')
#     cs<-bdag.to.cdag(bs)
#     print(cs)
#     cat('\n')
#   }
#   #the equivalent function
#   cat('Equivalent pairs:\n');
#   for (j in 1:25 ) {
#     for (i in 1:25) { 
#       if ( equivalent(cdag.to.bdag(A[j,]),cdag.to.bdag(A[i,])) ) 
#         print(c(j,i))
#     }
#   }
#   cat('\n')
# 
#   cat('Equal pairs:\n');
#   for (j in 1:25 ) {
#     for (i in 1:25) { 
#       if ( equal(cdag.to.bdag(A[j,]),cdag.to.bdag(A[i,])) ) 
#         print(c(j,i))
#     }
#   }
#   cat('\n')
# 
#   cat('BDAG:\n')
#   print(B<-array(c(0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0),c(4,4)))
#   cat('has equivalence class:\n')
#   print(M<-bdag.to.mdag(B))
#   cat('with these members:\n')
#   print(bdag.to.cdag(mdag.to.bdag(M)))
#   cat('\n')
# 
#   cat('BDAG:\n')
#   print(B<-array(c(0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0),c(4,4)))
#   cat('has equivalence class:\n')
#   print(M<-bdag.to.mdag(B))
#   cat('with these members:\n')
#   print(bdag.to.cdag(mdag.to.bdag(M)))
#   cat('\n')
# 
#   cat('BDAG:\n')
#   print(B<-array(c(0,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0),c(4,4)))
#   cat('has cycles:\n')
#   print(cycles(B))
#   cat('has v-structures:\n')
#   print(vstructures(B))
#   cat('has equivalence class:\n')
#   print(M<-bdag.to.mdag(B))
#   cat('with these members:\n')
#   print(bdag.to.cdag(mdag.to.bdag(M)))
#   cat('\n')
# 
#   cat('BDAG:\n')
#   print(B<-array(c(0,1,1,0,0,0,1,1,0,0,0,1,0,0,0,0),c(4,4)))
#   cat('has cycles:\n')
#   print(cycles(B))
#   cat('has v-structures:\n')
#   print(vstructures(B))
#   cat('has equivalence class:\n')
#   print(M<-bdag.to.mdag(B))
#   cat('with these members:\n')
#   print(bdag.to.cdag(mdag.to.bdag(M)))
#   cat('\n')
# 
#   cat('BDAG:\n')
#   print(B<-array(c(0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,0),c(4,4)))
#   cat('has cycles:\n')
#   print(cycles(B))
#   cat('has v-structures:\n')
#   print(vstructures(B))
#   cat('has equivalence class:\n')
#   print(M<-bdag.to.mdag(B))
#   cat('with these members:\n')
#   print(mdag.to.bdag(M))
#   print(bdag.to.cdag(mdag.to.bdag(M)))
#   cat('\n')
cat('testing orientation rule 1:\n');
B<-array(c(0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0),c(4,4))
cat('BDAG:\n');
print(B)
M<-bdag.to.mdag(B)
cat('MDAG:\n');
print(M)
cat('expanded:\n');
print(mdag.to.bdag(M))

cat('testing orientation rule 2:\n');
B<-array(c(0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0),c(5,5))
cat('BDAG:\n');
print(B)
cat('vstructures:\n');
print(vstructures(B))


M<-bdag.to.mdag(B)

cat('MDAG:\n');
print(M)
cat('expanded:\n');
print(mdag.to.bdag(M))
}