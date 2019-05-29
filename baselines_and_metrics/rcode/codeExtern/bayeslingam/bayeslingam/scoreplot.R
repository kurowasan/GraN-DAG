#for plotting the score calibration plot
#under the function used for the figure in the paper 
#this one 


#use for example pdf(file='scoreplot.pdf',width=16.5,height=1/2*9/13*16.5)

scoreplot<-function(dir='plotdata',by=1) {
  #K<-(dataindex-1)
  K<- length(dir(sprintf('%s/',dir),pattern="p[[:digit:]]*.Rdata"))#how many data sets
  Ks<-seq(from=1,to=K,by=by)

  bins<-c(5,15,25,35,45,55,65,75,85,95)

  abins<-array(0,c(2,length(bins)))
  bbins<-array(0,c(2,length(bins)))
  gbins<-array(0,c(2,length(bins)))
  lbins<-array(0,c(2,length(bins)))
  cbins<-array(0,c(2,length(bins)))

  #with a prior
 # adjust<-array(1,c(101,2))

  for ( fileindex in 1:length(Ks) ) {
    #cat('index:');print(i);
    #load parameters
    load(sprintf('%s/p%i.Rdata',dir,Ks[fileindex]));

    load(sprintf('%s/a%i.Rdata',dir,Ks[fileindex]));

    for ( i in 1:nrow(results$DAGs) ) { #this corrects dagind to a correct index
        if ( equal( parameters$DAG, cdag.to.bdag(results$DAGs[i,])  ) ) {
          parameters$dagind<-i
          break;
        }
     }

    i<-which.min(abs(bins-100*results$prob[parameters$dagind]))
    abins[1,i]<-abins[1,i]+1
    for (p in results$prob[-parameters$dagind]) {
      i<-which.min(abs(bins-100*p))
      abins[2,i]<-abins[2,i]+1
    }


    load(sprintf('%s/b%i.Rdata',dir,Ks[fileindex]));

    for ( i in 1:nrow(results$DAGs) ) { #this corrects dagind to a correct index
        if ( equal( parameters$DAG, cdag.to.bdag(results$DAGs[i,])  ) ) {
          parameters$dagind<-i
          break;
        }
     }

    i<-which.min(abs(bins-100*results$prob[parameters$dagind]))
    bbins[1,i]<-bbins[1,i]+1
    for (p in results$prob[-parameters$dagind]) {
      i<-which.min(abs(bins-100*p))
      bbins[2,i]<-bbins[2,i]+1
    }

    load(sprintf('%s/g%i.Rdata',dir,Ks[fileindex]));
    for ( i in 1:nrow(results$DAGs) ) { #this corrects dagind to a correct index
        if ( equal( parameters$DAG, cdag.to.bdag(results$DAGs[i,])  ) ) {
          parameters$dagind<-i
          break;
        }
     }


    i<-which.min(abs(bins-100*results$prob[parameters$dagind]))
    gbins[1,i]<-gbins[1,i]+1
    for (p in results$prob[-parameters$dagind]) {
      i<-which.min(abs(bins-100*p))
      gbins[2,i]<-gbins[2,i]+1
    }

    if ( any( fileindex == round(length(Ks)*c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))) ) {
      cat('.');
    }
  }
  cat('\n')

  predicted<-bins
  
  par(mfcol=c(1,3),oma=c(3,3,3,0))

  par(mfg=c(1,1))
  par(mai=c(0.5,0.5,0,0))
  happened<-100*abins[1,]/(abins[1,]+abins[2,])
  barplot(happened,ylim=c(0,100),
          names.arg=c("-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-"),
                      cex=1.5,cex.axis=1.8,cex.names=1.8)
 # title('BL-GL',cex=2.0)

  par(mfg=c(1,2))
   par(mai=c(0.5,0.5,0,0))
  happened<-100*bbins[1,]/(bbins[1,]+bbins[2,])
  barplot(happened,ylim=c(0,100),
          names.arg=c("-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-"),
                      cex=1.5,cex.axis=1.8,cex.names=1.8)
  #title('BL-MoG',cex=2.0)

  par(mfg=c(1,3))
  par(mai=c(0.5,0.5,0,0))
  happened<-100*gbins[1,]/(gbins[1,]+gbins[2,])
  barplot(happened,ylim=c(0,100),
          names.arg=c("-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-"),
          cex=1.5,cex.axis=1.8,cex.names=1.8)


   mtext("BL-GL",side=3,at=1/6,outer=TRUE,cex=1.5,line=1)
   mtext("BL-MoG",side=3,at=3/6,outer=TRUE,cex=1.5,line=1)
   mtext("GH",side=3,at=5/6,outer=TRUE,cex=1.5,line=1)

   mtext("% of true graphs",side=2,at=1/2,outer=TRUE,cex=1.5,line=1)
   mtext("predicted",side=1,at=1/2,outer=TRUE,cex=1.5,line=1)

}
#use pdf(file='scoreplot2.pdf',width=16.5,height=1/3*16.5)