#datadir<<-'/home/ajhyttin/svn/bayeslingam/trunk/personal/antti'

#this function loads all directories needed for the bayeslingam article project
loadbayeslingam <- function( ... ) {
  #if cauzality_part is foun
  if ( any( ls(envir=globalenv()) == "cauzality_path" ) ) {

    common_dir<-sprintf("%s/trunk/common",cauzality_path)
    bayeslingam_dir<-sprintf("%s/trunk/bayeslingam",cauzality_path)
    rbayeslingam_dir<-sprintf("%s/trunk/bayeslingam/rbayeslingam",cauzality_path)
    rlingam_dir<-sprintf("%s/trunk/rlingam",cauzality_path)
    rdeal_dir<-sprintf("%s/trunk/rdeal",cauzality_path)
    rpc_dir<-sprintf("%s/trunk/rpc",cauzality_path)
    etudag_dir<-sprintf("%s/trunk/ETUDAG",cauzality_path)

  }

  source(sprintf("%s/sourcedir.R",common_dir))
  sourceDir(bayeslingam_dir, trace=FALSE) #bayeslingam directory
  sourceDir(rbayeslingam_dir, trace=FALSE) #the algorithm directory
  sourceDir(rlingam_dir, trace=FALSE)
  sourceDir(rdeal_dir, trace=FALSE)
  sourceDir(rpc_dir, trace=FALSE)
  sourceDir(common_dir, trace=FALSE)
  sourceDir(etudag_dir, trace=FALSE)

}