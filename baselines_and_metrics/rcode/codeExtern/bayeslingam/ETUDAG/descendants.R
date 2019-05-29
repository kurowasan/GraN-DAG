# gives back the descendants of a set Z in a BDAG B
# or integer(0) if Z doesn't have descendents

descendants <- function(B, Z) {

  accB <- accessible(B)
  desc <- which(accB[,Z] != 0)%%ncol(B)
  if (0 %in% desc) desc <- sort(setdiff(union(ncol(B),desc),0))

  desc

}