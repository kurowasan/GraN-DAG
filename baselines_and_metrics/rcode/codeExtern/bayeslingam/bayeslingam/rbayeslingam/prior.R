#forwards the prior to the prior of the correct model

logprior<-function(p,
                   mixtures=2,
                   model='MoG',
                   derivatives=TRUE) {
  if ( model == 'MoG' ) {
    R<-logpriorMoG(p,mixtures=mixtures,derivatives=derivatives)
  } else if ( model =='GL') {
    R<-logpriorGL(p,derivatives=derivatives)
  }
  R
}



