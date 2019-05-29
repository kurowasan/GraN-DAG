# testbayeslingam - test and evaluate the bayeslingam method
#
# SYNTAX:
# res <- testbayeslingam( task=1, quicker=FALSE )
#
# INPUTS:
# task     - task number (which test to perform)
# quicker  - boolean specifying whether to do quick version of tests
#            (note that FALSE, the default, plots everything into files)
#
# OUTPUTS:
# (varies with the task)
#

testbayeslingam <- function( task=1 ) {

  # Show user we are embarking on testbayeslingam
  cat(sprintf('Testbayeslingam: task %i.\n',task))

  # This will hold any results output by the tests
  R<-NULL
  
  # Load all the required code
  loadbayeslingam()

  #-----------------------------------------------------------------------
  # SOME SIMPLE TESTS
  #-----------------------------------------------------------------------
  
  # --- Two observed variables, non-gaussian data ------------------------
  if ( task == 1 ) {  
    cat('Quick test of bayeslingam with MoG model, non gaussian data, 2 nodes.\n')
    R<-quicktestbayeslingam(drawgraph=TRUE)
  }

  # --- Three observed variables, non-gaussian data ----------------------
  if ( task == 2 ) {
    cat('Quick test of bayeslingam with MoG model, non gaussian data, 3 nodes.\n')
    R<-quicktestbayeslingam(drawgraph=TRUE,nvars=3,model='GL')
  }

  # --- Two observed variables, gaussian data ----------------------------
  if ( task == 3 ) {  
    cat('Quick test of bayeslingam with GL model, gaussian data 2 nodes.\n')
    R<-quicktestbayeslingam(drawgraph=TRUE,nvars=2,nongaussiandata=FALSE,model='GL')
  }

  # --- Three observed variables, gaussian data --------------------------
  if ( task == 4 ) {
    cat('Quick test of bayeslingam with MoG, gaussian data 3 nodes.\n')
    R<-quicktestbayeslingam(drawgraph=TRUE,nvars=3,nongaussiandata=FALSE)
  }

  # --- Five variables greedy search
  if ( task == 5 ) {
    cat('Testing greedy search with 5 variables.\n')
    cat('correctDAG:\n');print(D$parameters$DAG)
    D<-createCaseData(list(nvars=5,N=5000))
    R<-greedybayeslingam(D$X)
    cleanResults(R)
    cat('Extending the result a bit.\n')
    R2<-extendGreedyResult(R,D$X)
    cleanResults(R)
  }

  if ( task == 6 ) {
    cat('Testing bayeslingam with MCMC.\n')
    D<-createData(list(nvars=2,N=5000,logk=0.3))
    cat('correctDAG:\n');print(D$parameters$DAG)

    R<-bayeslingam(D$X,model='GL',mcmc=TRUE)
    cleanResults(R)
  }

  #-----------------------------------------------------------------------
  # THOROUGH TESTS WITH VARIOUS NON-GAUSSIANITIES AND SAMPLE SIZES
  #-----------------------------------------------------------------------
  
  # --- Generate simulated data for 2 variables --------------------------
  
  if ( task == 7 ) {

    # Tell the user what we are doing
    cat('Creating 117000 data set for the square plot',
        'and score plot with 2 node networks.\n')
    wait(10)

    # For reproducibility, set the seed to zero
    set.seed(0)

    # Logarithm of exponent for non-gaussianity
    logks <- seq(from=-1,to=1,by=0.25)

    # Runs to average over
    runs <- 1:1000

    # Sample sizes to try
    samples <- round(10^seq(from=1,to=4,by=0.25))

    # This actually creates the dataset
    createDataSet( dir='plotdata', logks=logks, runs=runs,
                   samples=samples, nvars=2 )

  # --- Run all algorithms on the two-variable datasets ------------------
    cat('Running algorithms for the data set of 2 node networks.\n')
    cat('WARNING: THIS TASK WILL TAKE SEVERAL DAYS.\n')
    wait(10)
    runtestbayeslingam(c(0,1,2,3,4),dir='plotdata')


    
  }

  # --- Generate simulated data for 3 variables --------------------------
  
  if ( task == 8 ) {

    # Tell the user what we are doing
    cat('Creating 117000 data set for the square plot',
        'and score plot with 3 node networks.\n')
    wait(10)
    
    # For reproducibility, set the seed to zero
    set.seed(0)

    # Logarithm of exponent for non-gaussianity    
    logks<-seq(from=-1,to=1,by=0.25)
    
    # Runs to average over
    runs<-1:1000

    # Sample sizes to try    
    samples<-round(10^seq(from=1,to=4,by=0.25))

    # This actually creates the dataset    
    createDataSet( dir='plotdata3', logks=logks, runs=runs,
                   samples=samples, nvars=3 )

    cat('Running algorithms for the data sets of 3 node networks.\n')
    cat('LONG RUN SELECTED. 117000 DATA SETS.\n')
    cat('WARNING: THIS TASK WILL TAKE SEVERAL DAYS.\n')
    wait(10)
    runtestbayeslingam(c(0,1,2,3,4),dir='plotdata3')

  }

  # --- Plotting the 'square' figure for 2 variables  --------------------
  
  if ( task == 9 ) {
    cat('Plotting the square figure for 2 node networks.\n')
    wait(5)
    pdf(file='squareplot2.pdf',width=16.5,height=9/13*16.5)
    squareplot(dir='plotdata')
    dev.off()
    cat('See squareplot2.pdf\n')
  }

  # --- Plotting the 'square' figure for 3 variables  --------------------
  
  if ( task == 10 ) {
    cat('Plotting the square figure for 3 node networks.\n')
    wait(5)
    pdf(file='squareplot3.pdf',width=16.5,height=9/13*16.5)
    squareplot(dir='plotdata3')
    dev.off()
    cat('See squareplot3.pdf\n')
  }

  # --- Calibration plots for 2 variables (Bayeslingam only!)  -----------

  if ( task == 11 ) {
    cat('Plotting the scoreplot-figure, only Bayeslingam, 2 node networks.\n')
    wait(5)
    pdf(file='scoreplot2.pdf',width=16.5,height=1/2*9/13*16.5)
    scoreplot(dir='plotdata',by=1)
    dev.off()
    cat('See scoreplot2.pdf\n')
  }

  #-----------------------------------------------------------------------
  # TESTS WITH PSEUDOREAL DATA 
  #-----------------------------------------------------------------------

  # --- Creating pseudoreal data -----------------------------------------
  
  if ( task == 12 ) {
    cat('Creating pseudoreal datasets.\n')
    wait(10)
    set.seed(0)
    createPseudoDataSets()

  # --- Running all algorithms on the pseudoreal data --------------------
    cat('Running algorithms on the pseudoreal datasets.\n')
    wait(10)
 
    for  (i in 1:20 ) {
      runtestbayeslingam(dir=sprintf('pseudodata%i',i),pseudo=TRUE)
    }    
  }

  # --- Plotting the results ---------------------------------------------

  if (task == 13 ) {
    cat('Plotting the pseudorealplot-figure, all algorithms, pseudoreal data.\n')
    wait(5)
    pdf(file='pseudorealplot.pdf',width=16.5,height=1/3*16.5)
    pseudoRealPlot()
    dev.off()
    cat('See pseudorealplot.pdf\n')
  }

  if ( task == 14 ) {
    # Tell the user what we are doing
    cat('Creating 1170 data set for the square plot',
        'and score plot with 4 node networks.\n')
    wait(10)
    
    # For reproducibility, set the seed to zero
    set.seed(0)

    # Logarithm of exponent for non-gaussianity    
    logks<-seq(from=-0.5,to=0.5,by=0.25)
    
    # Runs to average over
    runs<-1:10

    # Sample sizes to try    
    samples<-round(10^seq(from=1,to=3,by=0.25))

    # This actually creates the dataset    
    createDataSet( dir='plotdata4', logks=logks, runs=runs,
                   samples=samples, nvars=4 )

    cat('Running algorithms for the data sets of 4 node networks.\n')
    cat('WARNING: THIS TASK WILL TAKE SEVERAL HOURS.\n')
    wait(10)
    #skipping GH-dealer, the use of it is not ok for this many nodes
    runtestbayeslingam(c(0,1,3,4),dir='plotdata4')

    
  }

  if ( task == 15 ) {
    # Tell the user what we are doing
    cat('Creating 1170 data set for the square plot',
        'and score plot with 5 node networks.\n')
    wait(10)
    
    # For reproducibility, set the seed to zero
    set.seed(0)

    # Logarithm of exponent for non-gaussianity    
    logks<-seq(from=-0.5,to=0.5,by=0.25)
    
    # Runs to average over
    runs<-1:10

    # Sample sizes to try    
    samples<-round(10^seq(from=1,to=3,by=0.25))

    # This actually creates the dataset    
    createDataSet( dir='plotdata5', logks=logks, runs=runs,
                   samples=samples, nvars=5 )

    cat('Running algorithms for the data sets of 5 node networks.\n')
    cat('WARNING: THIS TASK WILL TAKE SEVERAL HOURS.\n')
    wait(10)
        #skipping GH-dealer, the use of it is not proper for this many nodes
    runtestbayeslingam(c(0,1,3,4),dir='plotdata5')

  }

  if ( task == 16 ) {
    cat('Plotting the squareplot for 4 variables.\n')
    squareplot2('plotdata4')
  }

  if ( task == 17 ) {
    cat('Plotting the squareplot for 4 variables.\n')
    squareplot2('plotdata4')
  }

  if (task == 18 ) { #this produces the case
    caseplot()
  }

  # Do not return any results
  #R
}
