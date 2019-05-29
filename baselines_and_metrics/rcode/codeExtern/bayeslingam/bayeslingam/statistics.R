statistics<-function(dir='plotdata',by=1,interval=FALSE) {
  nvars<-3
  algos<-5

  logk<-seq(from=-1,to=1,by=0.25)  #this should be in y-axis
  samples<-seq(from=1,to=4,by=0.25) #this should be in y-axis


  # the number of datasets
  K<- length(dir(sprintf('%s/',dir),pattern="a[[:digit:]]*.Rdata"))

  cat('Total of testruns used:');print(K/by)

  Ks<-seq(from=1,to=K,by=by)

  logscores<-array(NA,c(algos,length(Ks)))
  quadscores<-array(NA,c(algos,length(Ks)))
  binscores<-array(NA,c(algos,length(Ks)))
  classscores<-array(NA,c(algos,length(Ks)))

  for ( fileindex in 1:length(Ks) ) {
    cat('index:');print(fileindex);
    #load parameters
    load(sprintf('%s/p%i.Rdata',dir,Ks[fileindex]));

    letter <- c( "a","b","g","l","c" )

    for ( k in 1:algos) {
        load(sprintf('%s/%s%i.Rdata',dir,letter[k],Ks[fileindex]));

        rating<-rate(results,correctDAG=parameters$DAG)

        logscores[k,fileindex]<-rating$log
        quadscores[k,fileindex]<-rating$quad
        binscores[k,fileindex]<-rating$binary
        classscores[k,fileindex]<-rating$class
    }

    if ( any( fileindex == round(length(Ks)*
        c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))) ) {
      cat('.');
    }
  }
  cat('\n');
  cat('bincores:\n');
  print(t(binscores));

  means<-rbind(apply(logscores,1,mean),apply(quadscores,1,mean),
               apply(binscores,1,mean),apply(classscores,1,mean) )

  if ( interval ) {
    sdd<-array(NA,c(nrow(means),ncol(means)))

    sdd[1,1:5]<-apply(logscores[1:5,],1,sd,na.rm=TRUE)
    sdd[2,1:5]<-apply(quadscores[1:5,],1,sd,na.rm=TRUE)
    sdd[3,1:5]<-apply(binscores[1:5,],1,sd,na.rm=TRUE)
    sdd[4,1:5]<-apply(classscores[1:5,],1,sd,na.rm=TRUE)
    sdd<-sdd/sqrt(ncol(logscores))
    sdpos<-means+sdd
    sdneg<-means-sdd

  } else {
    sdpos<-NA
    sdneg<-NA
  }

   algo<-c("BL-GL","BL-MoG","GH","Lingam","PC")
   score<-c("log","quad","bin","class" )

   D<-data.frame(cbind(score,means) )
   names(D)<-c("score",algo)
   print(D)
   print(sdd)
 # quadscores
   list(means=means,sdpos=sdpos,sdneg=sdneg)
}

#use pdf(width=16.5,height=4*16.5)   
#its inches but a good picture when just give the centimeters

