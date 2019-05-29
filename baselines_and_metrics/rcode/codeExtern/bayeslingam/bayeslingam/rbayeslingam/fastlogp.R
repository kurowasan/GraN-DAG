# fastlogp:
# calculates the score of any single DAG from the precomputed component scores 
# the component scores have to be calculated
#
# SYNTAX:
# score<-fastlogp( DAG, components )
#
# INPUT:
# DAG          - DAG vector
# components   - precomputed scores for the components
#
# OUTPUT:
# score        - score for the dag, added up from the component scores

fastlogp <- function( DAG, components )
{

  # transform the dag to a B matrix
  Bup <- cdag.to.bdag(DAG)

  # initialize the score variable
  score<-0

  # go through the rows of B that correspond to one node and its parents
  for (i in (1:nrow(Bup)) ) {

    # find the appropriate component and add its score to the DAG score
    for (j in (1:length(components$node)) ) {
      if ( components$node[j] == i && #correct node
           all( components$edges[j,] == Bup[i,] ) ) { #correct parent conf. 
        score <- score+components$score[j]
        break;  #once found move to the next row
      }
    }
  }
  # return the DAG score
  score
}
