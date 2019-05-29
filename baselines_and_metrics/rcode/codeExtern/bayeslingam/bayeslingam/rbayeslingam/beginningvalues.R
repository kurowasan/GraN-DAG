# beginningvalues - calculates initial values for numerical optimization
#
# SYNTAX:
# beginningvalues( y, X, mixtures, model )
#
# INPUTS:
# y        - the values of the child
# X        - the values of the parents (must be matrix N x number of parents)
# mixtures - the number of components for the 'mixnorm' model
# model    - 'MoG' or 'GL'
#
# OUTPUT:
# p        - vector of starting parameter values
#


beginningvalues<-function(y,X,mixtures,model) {
  if ( model == 'MoG') {
    R<-beginningvaluesMoG(y,X,mixtures)
  } else if ( model == 'GL' ) {
  }

  R
}