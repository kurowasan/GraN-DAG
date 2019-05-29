#SHOULD THIS BE A PATTERN???

#are the constructors really necessary? well, for the equal operator they are

#mdag<-function( M ) { #makes an mdag from a matrix M
#  attr(M,'class')<-'MDAG'
#  mode(M)<-'integer' #for storing!
#  M
#}

# THIS IS THE DEFINITION OF A MDAG FROM THE PC ALGORITHM PACKAGE
#within MDAG matrix
#MDAG[i,j]=0 if there is no edge between the variables
#MDAG[i,j]=2 if there is edge i->j
#MDAG[i,j]=1 if there is edge i<-j
#MDAG[i,j]=3 if there is edge i-j