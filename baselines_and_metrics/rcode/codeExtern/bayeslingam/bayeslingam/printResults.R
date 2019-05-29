printResults<-function() {
  algos<-4
  scores<-4

  dags<-alldags(nvars)

  K<- length(dir(sprintf('%s/',dir),pattern="p[[:digit:]]*.Rdata"))#how many data sets
  cat('Total of testruns used:');print(K/by)

  Ks<-seq(from=1,to=K,by=by)

  N<-0
  scores<-array(0,c(algos,scores))

  for ( fileindex in 1:length(Ks) ) {
    #cat('index:');print(fileindex);
    #load parameters
    N<-N+1
    load(sprintf('%s/p%i.Rdata',dir,Ks[fileindex]));

    letter <- c( "b","g","l","c" )
    for ( k in 1:algos) {
      load(sprintf('%s/%s%i.Rdata',dir,letter[k],Ks[fileindex]));
      rating<-rate(results,correctIndex=parameters$dagind,dags=dags)
      scores[k,1]<-scores[k,1]+rating$log
      scores[k,2]<-scores[k,2]+rating$quad
      scores[k,3]<-scores[k,3]+rating$binary
      scores[k,4]<-scores[k,4]+rating$class
    }

    if ( any( fileindex == round(length(Ks)*c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))) ) {
      cat('.');
    }
  }
  cat('\n')

  cat('Check that all dags are equally likely:\n')
  cat('N:');print(N)
  cat('dagind:');print(D)

  #now we would like to normalize each score to 0:1 in such a way that the plot shades are equivalent
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

  algo<-c("BL","Deal","Lingam","PC")
  score<-c("log","quad","bin","class" )
  for (i in 1:4) { # i is the score
    for (j in 1:algos) { #j is the algorithm
      if ( i == 1 ) {
        Z<-logscores[,,j]
      } else if (i == 2 ) {
        Z<-quadscores[,,j]
      } else if (i == 3) {
        Z<-binscores[,,j]
      } else if (i == 4 ) {
        Z<-classscores[,,j]
      }
       par(mfg=c(i,j))

      vmar<-rep(0.0,4)
      lx<-'n'
      ly<-'n'
      xl<-''
      yl<-''

       if ( i == 4 ) {
         #lx<-NULL
         #xl<-'log10(samples)'
         #vmar[1]<-4
       }
       if ( j == 1 ) {
         #ly<-NULL
         #yl<-'log k'
         #vmar[2]<-4
       }
      par(mar=vmar,bty="n")

      image(x=samples,y=logk,z =Z, col=gray((0:32)/32),zlim=c(0,1),xlab=xl,ylab=yl,xaxt=lx,yaxt=ly,bty="o")
       if ( i == 4 ) {
         axis(side=1,outer=TRUE,labels=c("10"," ","100"," ","1000"," ","10000"),at=c(1,1.5,2,2.5,3,3.5,4))
         #lx<-NULL
         #xl<-'log10(samples)'
         #vmar[1]<-4
       }
       if ( j == 4 ) {
         axis(side=4,outer=TRUE)
         #ly<-NULL
         #yl<-'log k'
         #vmar[2]<-4
       }

      #title(sprintf('%s:%s',algo[j],score[i]));
    }
  }
 # mtext("Algorithms:Bayeslingam,Deal,Lingam,PC",at=c(   ),outer=TRUE)

   mtext("Bayeslingam",side=3,at=1/8,outer=TRUE,cex=1.3,line=1)
   mtext("Deal",side=3,at=3/8,outer=TRUE,cex=1.3,line=1)
   mtext("Lingam",side=3,at=5/8,outer=TRUE,cex=1.3,line=1)
   mtext("PC",side=3,at=7/8,outer=TRUE,cex=1.3,line=1)

   mtext("log",side=2,at=7/8,outer=TRUE,cex=1.3,line=1)
   mtext("quadratic",side=2,at=5/8,outer=TRUE,cex=1.3,line=1)
   mtext("binary",side=2,at=3/8,outer=TRUE,cex=1.3,line=1)
   mtext("class",side=2,at=1/8,outer=TRUE,cex=1.3,line=1)

  mtext("samples",side=1,at=1/2,outer=TRUE,cex=1.3,line=3)
  mtext("log k",side=4,at=1/2,outer=TRUE,cex=1.3,line=3)

}