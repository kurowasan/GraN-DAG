# subfunction of dseparations
# returns TRUE if the x and y are d-connected given Z
# x != y, x,y not in Z (x,y and all elements in Z are integers)
# B = BDAG matrix
# recursive function, income denotes if the arrow from the previous step is
# pointing into x (TRUE) or not.
# var is the set of all variables already visited (to not get in endless loop)

dconnected <- function(B, x, y, Z, income, fromx, var=NULL) {

  if (x == y) return(TRUE)

  dcon <- FALSE

  if (income) {
    
    # have fromx -> x
    
    # connecting path fromx -> x <- parent(x) if x in Z or desc(x) in Z OR ...
    if (x%in%Z | any(descendants(B,x)%in%Z)) { # open path through collider
      x.par <- setdiff(which(B[x,]!=0),fromx)
      # take fromx out of parents, came from fromx, don't want to go back there
      for (j in x.par) {
        if (!j%in%var) dcon <- dconnected(B,j,y,Z,FALSE,x,union(var,fromx))
        if (dcon) return(TRUE)
      }
    } # end if x in Z
    
    # ... OR connecting path fromx -> x -> child(x) if x not in Z
    if (!x%in%Z) { # open path through chains
      x.chi <- which(B[,x]!=0)
      for (j in x.chi) {
        if (!j%in%var) dcon <- dconnected(B,j,y,Z,TRUE,x,union(var,fromx))
        if (dcon) return(TRUE)
      }
    } # end if x in Z
  
  } # end income-if

  else { # income = FALSE
    
    # have fromx <- x
    
    # no connecting path fromx <- x <- parent(x) if x in Z and
    # no connecting path fromx <- x -> child(x) if x in Z OR ...
    if (x%in%Z) return(FALSE)

    # ... OR connecting path fromx <- x <- parent(x) if x in Z and
    # also connecting path fromx <- x -> child(x) if x in Z
    else { # x not in Z, open path through chain or fork
      x.par <- which(B[x,]!=0)
      x.chi <- setdiff(which(B[,x]!=0),fromx)
      for (j in x.par) {
        if (!j%in%var) dcon <- dconnected(B,j,y,Z,FALSE,x,union(var,fromx))
        if (dcon) return(TRUE)
      }
      for (j in x.chi) {
        if (!j%in%var) dcon <- dconnected(B,j,y,Z,TRUE,x,union(var,fromx))
        if (dcon) return(TRUE)
      }
    } # end else x not in Z

  } # end income-else

  if (!dcon) return(FALSE)

}