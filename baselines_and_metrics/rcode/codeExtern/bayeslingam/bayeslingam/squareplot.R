#Plot the squareplot of the paper
#with by=10 or 100 a faster result
#saves files for the results

squareplot<-function(dir='plotdata',by=1,calc=TRUE) {

  algos<-5
 # require(grDevices) # for colours PROBABLY NOT NEEDED
  logk<-seq(from=-1,to=1,by=0.25)  #this should be in y-axis
  samples<-seq(from=1,to=4,by=0.25) #this should be in y-axis

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

  if ( calc ) {
    for ( fileindex in 1:length(Ks) ) {
      #cat('index:');print(fileindex);
      #load parameters
      load(sprintf('%s/p%i.Rdata',dir,Ks[fileindex]));
      #print(c(parameters$logk,parameters$N))
      j<-which.min(abs(logk-parameters$logk))
      i<-which.min(abs(round(10^samples)-parameters$N))
      #print(c(i,j))
      N[i,j]<-N[i,j]+1

      letter <- c( "a","b","g","l","c" )

      for ( k in 1:algos) {
        filename=sprintf('%s/%s%i.Rdata',dir,letter[k],Ks[fileindex]);

        if ( file.access(filename,0) != 0 ) {#this means that the file is not found
          #do nothing!!!
          
        } else {
          load(filename) 
        }
        #lines added to fit the old to the newer result structure
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
  
    #this saves or loads these files that take quite a bit of calculation also
    #from the already saved results
    save(logscores,file='logscores.Rdata')
    save(binscores,file='binscores.Rdata')
    save(quadscores,file='quadscores.Rdata')
    save(classscores,file='classscores.Rdata')
  } else {
  #load these if everything is already calculated!
    load(file='logscores.Rdata')
    load(file='binscores.Rdata')
    load(file='quadscores.Rdata')
    load(file='classscores.Rdata')
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

  algo<-c("BL-GL","BL-MoG","G-H","LiNGAM","PC")
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

      image(x=samples,y=logk,z =Z,
            col=gray((0:32)/32),zlim=c(0,1),xlab=xl,ylab=yl,
            xaxt=lx,yaxt=ly,bty="o")
       if ( i == 4 ) {
          if ( j == 1 ) {
            axis(side=1,outer=TRUE,labels=c("10"," ","100"," ","1000",
                  " ","10000"),at=c(1,1.5,2,2.5,3,3.5,4),cex.axis=1.8)
          } else {
            axis(side=1,outer=TRUE,labels=c(""," ","100"," ","1000",
                 " ","10000"),at=c(1,1.5,2,2.5,3,3.5,4),cex.axis=1.8)
          }

       }
       if ( j == 4 ) {
         axis(side=4,outer=TRUE,labels=c("1/2","1","2"),
              at=c(log(1/2),log(1),log(2)),cex.axis=1.8)

       }


    }
  }

   mtext("BL-GL",side=3,at=1/10,outer=TRUE,cex=1.5,line=1)
   mtext("BL-MoG",side=3,at=3/10,outer=TRUE,cex=1.5,line=1)
   mtext("GH",side=3,at=5/10,outer=TRUE,cex=1.5,line=1)
   mtext("LiNGAM",side=3,at=7/10,outer=TRUE,cex=1.5,line=1)
   mtext("PC",side=3,at=9/10,outer=TRUE,cex=1.5,line=1)

   mtext("binary",side=2,at=7/8,outer=TRUE,cex=1.5,line=1)
   mtext("class",side=2,at=5/8,outer=TRUE,cex=1.5,line=1)
   mtext("log",side=2,at=3/8,outer=TRUE,cex=1.5,line=1)
   mtext("quadratic",side=2,at=1/8,outer=TRUE,cex=1.5,line=1)

  mtext("samples",side=1,at=1/2,outer=TRUE,cex=1.5,line=3)
  mtext("q",side=4,at=1/2,outer=TRUE,cex=1.5,line=3)

}

#use pdf(file='squareplot2.pdf',width=16.5,height=4/5*9/13*16.5) 
#  its inches but a good picture when just give the centimeters

