# dseparations returns all the d-separations in the DAG
# possibility to include hiddenvariables in DAG

# need package "gtools" for fct combinations

dseparations <- function(BDAG, HVmatrix=NULL, output_info=FALSE) {

  library(gtools)

  O <- BDAG # DAG over observed variables
  L <- HVmatrix # matrix of hidden/latent variables, or NULL
  nvar <- nrow(O)

  B <- bdaghidd.to.hdag(O,L) # DAG over observed and hidden variables, if no
                             # hidden variables, this is simply O

  # ### get sets to test: x indep y given Z ### #
  # x,y, Z in observed variables
  allpairs <- combinations(n=nvar,r=2) # points x and y

  # cut out all pairs of nodes which have an edge between each other because
  # these nodes can never be d-separated (so don't have to check)!!
  for (i in 1:nvar) {
    temp <- which(B[,i]==1)
    for (j in temp) {
      pair <- sort(c(i,j))
      ind <- intersect(which(allpairs[,1]==pair[1]),
                       which(allpairs[,2]==pair[2]))
      allpairs <- allpairs[-ind,]
      if (is.null(dim(allpairs))) allpairs <- t(as.matrix(allpairs))
    }
  }

  # ### get d-separations ### #
  # x d-sep y given Z in the DAG over observed AND hidden varibales
  # check d-separations with 0 elements in Z, then with 1 element, then 2...
  # if x, y d-separated given some Z, take x,y out of allpairs list

  dsep <- array(0,dim=c(nrow(allpairs)*(nvar-1),nvar))
  cnt <- 0

  for (i in 0:(nvar-2)) { # i = number of elements in Z

    if (nrow(allpairs) == 0) break
    allsets <- allZsets(allpairs,nvar,i) # passible sets Z to x,y with i elem.

    for (j in 1:nrow(allsets)) {

      temp <- dseparated( B, allsets[j,1], allsets[j,2],
                  allsets[j,index(3,length(allsets[j,]))] )

      # if d-separated then temp = TRUE:
      if (temp) {

        # take pair x,y out of allpairs -> don't check with more elements in Z
        pair <- allsets[j,c(1,2)]
        ind <- intersect(which(allpairs[,1]==pair[1]),
                         which(allpairs[,2]==pair[2]))
        allpairs <- allpairs[-ind,]
        if (is.null(dim(allpairs))) allpairs <- t(as.matrix(allpairs))

        # output and save result
        if (output_info) {
          cat(allsets[j,1], 'independent', allsets[j,2], 'given',
              allsets[j,index(3,length(allsets[j,]))], '\n')
        }
        cnt <- cnt + 1
        dsep[cnt,1:(2+length( allsets[j, index(3,length(allsets[j,]))] ))] <-
            c( allsets[j,1],allsets[j,2],
               allsets[j,index(3,length(allsets[j,]))] )

      }

    } # end j

  } # end i

  # organize dsep

  dsep <- dsep[index(1,cnt),]

  # if dsep has 1 row only then it is a vector: convert it to matrix and return
  if (is.null(dim(dsep))) return(t(as.matrix(dsep)))

  # if dsep is the empty matrix: return
  if (dim(dsep)[1] == 0) return(dsep)

  # if dsep has more than 1 row: sort it after the first two columns and return
  if (!is.null(dim(dsep))) {
    ord <- order(dsep[,2])
    dsep <- dsep[ord,]
    ord <- order(dsep[,1])
    dsep <- dsep[ord,]
    return(dsep)
  }

}


allZsets <- function(allpairs,nvar,ncond) {

  # ### get sets to test: x indep y given Z ### #
  # x,y, Z in observed variables
  # x,y are given in allpairs, get all possible sets Z with ncond elements

  if (ncond==0) { # no elements in Z
    return(allpairs)
  }

  else { # ncond!=0 elements in Z
    possZ <- combinations(n=nvar,r=ncond)
    allsets <- array(0,dim=c(nrow(possZ)*nrow(allpairs),2+ncond))
    cnt <- 0
    for (j in 1:nrow(allpairs)) {
      for (k in 1:nrow(possZ)) {
        if (length(intersect(possZ[k,],allpairs[j,])) == 0) {
          cnt <- cnt + 1
          allsets[cnt,] <- c(allpairs[j,],possZ[k,])
        } # end if
      } # end k
    } # end j
  } # end else

  if (is.null(dim(allsets[index(1,cnt),]))) {
    allsets <- as.matrix(allsets[index(1,cnt),])
    dim(allsets) <- c(1,2+ncond)
    return(allsets)
  }

  allsets[index(1,cnt),]

}