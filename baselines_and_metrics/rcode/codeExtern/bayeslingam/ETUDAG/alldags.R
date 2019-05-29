#---------------------------------------------------------------------------
# alldags:
# returns a list of all DAGs over n variables
# the DAGs are represented as follows:
# each DAG is one row of the returned matrix
# the first n elements of each row denotes a causal order of the given DAG
# the remaining elements give the elements of the corresponding lower
# triangular matrix B, where element bij=1 if variable j is a parent of
# variable i, _in the causal order_.
#
# examples might make this a bit clearer:
# say n=3. here are some examples:
# 2 3 1 1 0 1 represents x2 -> x3 -> x1
# 1 3 2 0 1 1 represents x1 -> x2 <- x3
# 1 2 3 0 0 0 represents x1 , x2 , x3
# 3 2 1 1 0 0 represents x3 -> x2 , x1
# where ',' denotes non-adjacency
#
#---------------------------------------------------------------------------

alldags <- function( n ) {
  if ( n == 5 & any( ls(envir=globalenv()) == "cauzality_path" ) ) {
    #cat('Loading alldags(5) from a file.\n')
    load(sprintf('%s/trunk/ETUDAG/dags5.Rdata',cauzality_path))
    D<-dags
  } else if ( n == 6 ) {
    D<-pdags(n)
  } else {
    # initialize everything using the empty DAG
    #D <- matrix(0,1,n+n*(n-1)/2)
    ndags<-ap(n)$N;
    Dindex=0;
    D <- matrix(0,ndags,n+n*(n-1)/2)
    B <- matrix(0,ndags,n*n)
    #S <- rep(0,ndags)
  
    Dindex<-Dindex+1;
    D[1,1:n] <- 1:n
    #B[1,] OK already!
    #S ok as well
  
    # go through all permutations
    P <- all.perm( n )
    perms <- nrow(P)
    for (i in 1:perms) {
  
      # go through all connections
      for (j in 1:(2^(n*(n-1)/2))) {
  
        # convert decimal to binary
        binvec <- dec.to.bin(j,n*(n-1)/2)
        
        # set connection matrix
        Bmat <- matrix(0,n,n)
        Bmat[lower.tri(Bmat,diag=FALSE)] <- binvec
        
        # permute to appropriate permutation
        ip <- iperm(P[i,])
        Bmat <- Bmat[ip,ip]
        #Bmat <- Bmat[P[i,],P[i,]]
        
        # check if already in B, if so, go to next iteration
        Bvec <- Bmat
        dim(Bvec) <- c(1,n*n)
        alreadyinB <- FALSE
        #s<-sum(Bvec)
        #I<-which(s==S[1:Dindex])
        #if ( !is.null(I) ) {
          for (jj in 1:Dindex) {
            if (all(Bvec==B[jj,])) {
              alreadyinB <- TRUE
              break;
            }
          }
          if (alreadyinB) {
            next
          }
    
  
        # add it as a new DAG
        Dindex<-Dindex+1
        D[Dindex,] <- c(P[i,],binvec)
        B[Dindex,] <- Bvec
        #S[Dindex]<- s 
      }
    }#else
  }
  D
  
}

ap<-function( p ) {
  R<-list()
  R$N<-0
  for (k in 1:p) {
    R$N<-R$N+apk(p,k)$N
  }
  R
}

apk<-function(p, k) {
  #cat(sprintf('apk(%i,%i)\n',p,k))
  R<-list()
  R$N<-0
  if (p-k > 0) {
    for (n in 1:(p-k) ) {
      #cat(sprintf('apk(%i,%i),N=%i\n',p,k,n))
      #print(c((2^k-1)^n))
      #print((2^(k*(p-n-k))))
      #print(choose(p,k))
      #print(apk(p-k,n))
      R$N<-R$N+((2^k-1)^n)*(2^(k*(p-n-k)))*choose(p,k)*apk(p-k,n)$N
    }
  } else { #otherwise there is exactly one DAG with p nodes and k=p nodes with indegree 0
    R$N<-1
    #R$D<-as.matrix(t(c(1:p,rep(0,p))))#the empty DAG
  }
  #cat(sprintf('apk(%i,%i)=%i\n',p,k,S))
  R
}


