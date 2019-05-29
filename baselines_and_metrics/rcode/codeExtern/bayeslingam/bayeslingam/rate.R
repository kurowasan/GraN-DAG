#rating functions
#


rate<-function(results,correctDAG=NULL,probability_resolution=1e-5) {
  #this function rates the performance of bayeslingam with different criterions
  #results$DAGs is a matrix with each row as one 
  #results$prob is the vector of probabilities for each dags
    #cat('rate:\n');
    #print(correctDAG)

    if ( is.null(dim(correctDAG)) || nrow(correctDAG) == 1 ) {
      #cat('cdag.to.bdag:\n');
      correctDAG<-cdag.to.bdag(correctDAG)
    }
    #cat('rate continuing\n');

    losses<-list()
 
    #first find the winning index for class and binary score
    maxprob=max(results$prob);

    #get _ALL_ DAGs with equal probability
    winning_index=which( abs(results$prob-maxprob) < probability_resolution );

    #draw just one of them by random
    if ( length( winning_index) > 1 ) {
      winning_index=sample(winning_index,1);
    }
    #cat('rate continuing2\n');
    #print(winning_index)
    #print(results$DAGs)

    winning_bdag=cdag.to.bdag(results$DAGs[winning_index,]);


    # cat('rate continuing3\n');
    #cat('Winning BDAG:\n');
    #print(winning_bdag)

    #cat('Correct BDAG:\n');
    #print(correctDAG)

    #binary score, using ETUDAG package equal for bdags
    if ( equal( winning_bdag, correctDAG ) ) {
      losses$binary<-0;
    } else {
      losses$binary<-1;
    }

    #binary score, using ETUDAG package observationally equivalent for bdags
    #checks for v-structures and skeleton
    #cat('here:')
    if ( equivalent( winning_bdag, correctDAG ) ) {
      losses$class<-0;
    } else {
      losses$class<-1;
    }


    losses$log<-Inf
    losses$quad<-0
    #then log and quadratic score
    #here we dont really have the correct index!!!

    for ( i in 1:length(results$prob) ) { #go through all dags mentioned in the
                                          #result
      if ( equal( correctDAG, cdag.to.bdag(results$DAGs[i,])  ) ) {
        #now i is the correct index of previous
        losses$quad <- losses$quad + (results$prob[i]-1)^2
        losses$log<-(-1)*log(results$prob[i])
      } else {
        #now i is one of the 
        losses$quad <- losses$quad + (results$prob[i])^2
      }
      #the DAGs not mentioned in DAGs are assumed to have zero prob
      #if the correct DAG is not mentioned log score is infinitely automatically
      #i think these should work in all cases
    }#for i
    losses$arcsturned<-0
    losses$arcsadded<-0
    losses$arcsdeleted<-0
    losses$vsadded<-0
    losses$vsdeleted<-0

    #Calculating some structural properties of the winning bdag
    for (i in 1:nrow(winning_bdag) ) {
      for (j in 1:ncol(winning_bdag) ) {

        if ( winning_bdag[i,j] == 1  ) {
          if (correctDAG[i,j] == 1 ) {
            #OK!
          } else {
            if ( correctDAG[j,i] == 1 ) {
              losses$arcsturned<-losses$arcsturned + 1 
            } else {
              losses$arcsadded<-losses$arcsadded + 1
           } 
         }
        } else if ( winning_bdag[i,j] == 0 ) {
          if ( correctDAG[i,j] == 0 ) {
            #OK!
          } else if ( correctDAG[i,j] == 1 && winning_bdag[j,i] == 0 ) {
            losses$arcsdeleted<-losses$arcsdeleted+1
          }
        }
      }
    }
    correctvs<-vstructures(correctDAG)
    winningvs<-vstructures(winning_bdag)
    #print(correctvs)
    #print(winningvs)


    for ( i in index(1,nrow(winningvs)) ) {
      if ( !any( correctvs[,1] == winningvs[i,1] &
                correctvs[,2] == winningvs[i,2] &
                correctvs[,3] == winningvs[i,3]) ) {
          losses$vsadded<-losses$vsadded+1
      }
    }
    for ( i in index(1,nrow(correctvs)) ) {
      if ( !any( correctvs[i,1] == winningvs[,1] &
                correctvs[i,2] == winningvs[,2] &
                correctvs[i,3] == winningvs[,3]) ) {
          losses$vsdeleted<-losses$vsdeleted+1
      }
    }

  losses
}


rateObsolete<-function(results,correctIndex,probability_resolution=1e-5,dags=NULL) {
  #this function rates the performance of bayeslingam with different criterions
  #results$DAGs is a matrix with each row as one 
  #results$prob is the vector of probabilities for each dags

  losses<-list()
 
     #first log loss
    losses$log<-(-1)*log(results$prob[correctIndex])
  
    #then quadratic los
    losses$quad<-sum( (results$prob[-correctIndex])^2 ) + (results$prob[correctIndex]-1)^2
  
    #binary loss
    losses$binary<-0 #what if two have the same probability????
    #if ( max(round(100*results$prob)) > round(100*results$prob[correctIndex]) ) {
    #  losses$binary = 1
    #}
  
    m<-max(results$prob)
    i<-which( abs(results$prob-m) < probability_resolution ) 
    #the index of the probabilities that are as big as the best one
  
    if (length(i) >= 2 ) { 
      #if there are many with the same prob, draw one of them at random
      i<-sample(i,1)
    }
    if ( i == correctIndex ) {
      losses$binary<-0
    } else {
      losses$binary<-1
    }
  
  
    #then equivalence class binary? you get a point for the right equivalence class
    if (is.null(dags) ) { #do this only if dags are given
      dags<-alldags( max(results$DAGs ) ) #creating all dags
    }
  
    if ( equivalentDAGs(dags[correctIndex,],dags[i,] ) ) {
      losses$class <- 0
    } else {
      losses$class <- 1
    }
  losses
}