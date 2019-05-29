#This can be used for faster 4 and 5 variable squareplots

squareplot2<-function(dir='plotdata4',by=1) {

  algos<-4
 # require(grDevices) # for colours PROBABLY NOT NEEDED
  logk<-seq(from=-0.5,to=0.5,by=0.25)  #this should be in y-axis
  samples<-seq(from=1,to=3,by=0.25) #this should be in y-axis

#  r<-outer(logk,samples)
#  dags<-alldags(nvars)

  K<- length(dir(sprintf('%s/',dir),pattern="a[[:digit:]]*.Rdata"))
  #how many data sets 
  cat('Total of testruns used:');print(K/by)

  Ks<-seq(from=1,to=K,by=by)


  logscores<-array(0,c(length(samples),length(logk),algos))
  quadscores<-array(0,c(length(samples),length(logk),algos))
  binscores<-array(0,c(length(samples),length(logk),algos))
  classscores<-array(0,c(length(samples),length(logk),algos))
  
  
  
  N<-array(0,c(length(samples),length(logk)))
  

    for ( fileindex in 1:length(Ks) ) {
      #cat('index:');print(fileindex);
      #load parameters
      load(sprintf('%s/p%i.Rdata',dir,Ks[fileindex]));
      #print(c(parameters$logk,parameters$N))
      j<-which.min(abs(logk-parameters$logk))
      i<-which.min(abs(round(10^samples)-parameters$N))
      #print(c(i,j))
      N[i,j]<-N[i,j]+1
  
      letter <- c( "a","b","l","c" )
  
      for ( k in 1:algos) {
        load(sprintf('%s/%s%i.Rdata',dir,letter[k],Ks[fileindex]));

        rating<-rate(results,correctDAG=parameters$DAG)
  
        logscores[i,j,k]<-logscores[i,j,k]+rating$log
        quadscores[i,j,k]<-quadscores[i,j,k]+rating$quad
        binscores[i,j,k]<-binscores[i,j,k]+rating$binary
        classscores[i,j,k]<-classscores[i,j,k]+rating$class
      }
  
      #printing the progress
      if ( any( fileindex == round(length(Ks)*c(0.1,0.2,0.3,0.4,
                0.5,0.6,0.7,0.8,0.9,1.0))) ) {
        cat('.');
      }
    }


  cat('\n')

  cat('Check that all dags are equally likely:\n')
  cat('N:');print(N)

  #now we would like to normalize each score to 0:1
  # in such a way that the plot shades are equivalent
  #so set minimum to 0 and maximum to 1

  ma<-max(logscores)# THE WORST SCORE
  mi<-min(logscores) # THE BEST SCORE
  ma<-(-log(0.97))*sum(N) #LET 10% for the right dag be the worst score
  cat('log (min,max):');print(c(mi,ma))
  logscores[logscores > ma]<-ma 
  logscores<-(logscores-mi)/(ma-mi)


  ma<-max(quadscores)# THE WORST SCORE
  mi<-min(quadscores)# THE BEST SCORE
  cat('quad (min,max):');print(c(mi,ma))
  quadscores<-(quadscores-mi)/(ma-mi)
  

  ma<-max(binscores)# THE WORST SCORE
  mi<-min(binscores)# THE BEST SCORE
  cat('bin (min,max):');print(c(mi,ma))
  binscores<-(binscores-mi)/(ma-mi)

  ma<-max(classscores)# THE WORST SCORE
  mi<-min(classscores)# THE BEST SCORE
  cat('class (min,max):');print(c(mi,ma))
  classscores<-(classscores-mi)/(ma-mi)

  logscores<-1-logscores
  binscores<-1-binscores
  quadscores<-1-quadscores
  classscores<-1-classscores

#  xlabel<-c("10","","","","100","","","","1000","","","","10000")
  par(mfcol=c(4,algos),oma=c(5,3,3,5))

  algo<-c("BL-GL","BL-MoG","LiNGAM","PC")
  score<-c("log","quad","bin","class" )

  for (i in 1:4) { # i is the score
    for (j in 1:algos) { #j is the algorithm
      if ( i == 1 ) {
        Z<-binscores[,,j]
      } else if (i == 2 ) {
        Z<-classscores[,,j]
      } else if (i == 3) {
        Z<-logscores[,,j]
      } else if (i == 4 ) {
        Z<-quadscores[,,j]
      }
       par(mfg=c(i,j))

      vmar<-rep(0.0,4)
      lx<-'n'
      ly<-'n'
      xl<-''
      yl<-''

      par(mar=vmar,bty="n")

      image(x=samples,y=logk,z =Z,
            col=gray((0:32)/32),zlim=c(0,1),xlab=xl,ylab=yl,
            xaxt=lx,yaxt=ly,bty="o")

    }
  }

  mtext("binary",side=2,at=7/8,outer=TRUE,cex=1.5,line=1)
  mtext("class",side=2,at=5/8,outer=TRUE,cex=1.5,line=1)
  mtext("log",side=2,at=3/8,outer=TRUE,cex=1.5,line=1)
  mtext("quadratic",side=2,at=1/8,outer=TRUE,cex=1.5,line=1)

  mtext("samples",side=1,at=1/2,outer=TRUE,cex=1.5,line=3)
  mtext("q",side=4,at=1/2,outer=TRUE,cex=1.5,line=3)

}