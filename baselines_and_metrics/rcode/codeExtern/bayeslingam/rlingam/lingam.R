# lingam - perform the LiNGAM analysis
#
# SYNTAX:
# res <- lingam( X )
#
# INPUT:
# X     - Data matrix: each row is an observed variable, each
#         column one observed sample. The number of columns should
#         be far larger than the number of rows.
#
# OUTPUT:
# res$B     - Matrix of estimated connection strenghts
# res$stde  - Standard deviations of disturbance variables
# res$ci    - Constants
# res$k     - An estimated causal ordering
#
# (Note that B, stde, and ci are all ordered according to the same
#  variable ordering as X.)
#
# Version: 0.1 (8 Jan 2008) - first version for R
#
# See also ESTIMATE, PRUNE.

lingam <- function( X, output = FALSE ) {

  temp <- estimate( X, output )
  res <- prune( X, temp$k, output )
  res
  
}
