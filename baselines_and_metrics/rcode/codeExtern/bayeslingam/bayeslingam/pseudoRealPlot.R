pseudoRealGenerationPlot<-function( datasetindex=9 ) {
  X<-getRealData(datasetindex)
  parameters<-list(datasetindex=datasetindex)
  D<-createPseudoData(parameters)
  

  #x11()
  pdf(file='orginal9.pdf',width=16.5,height=9/13*16.5)
  pairs(X)
  dev.off()

  Y<-data.frame(D$X)
  names(Y)<-names(X)
  pdf(file='pseudo9.pdf',width=16.5,height=9/13*16.5)
  pairs(Y)
  dev.off()
}

pseudoRealPlot<-function(datadir='', datasetindex=c(1,3,4,6,9,10,13,16,19,20) ) {
  H<-array(NA,c(4,5,length(datasetindex)) ) #score, algorithm, datset
  Ineg<-array(NA,c(4,5,length(datasetindex)) ) #score, algorithm, datset
  Ipos<-array(NA,c(4,5,length(datasetindex)) ) #score, algorithm, datset

  for ( d in 1:length(datasetindex) ) {
    S<-statistics(sprintf('%s/pseudodata%i',datadir,datasetindex[d]),interval=TRUE,by=1)
    H[,,d]<-S$means
    Ipos[,,d]<-S$sdpos
    Ineg[,,d]<-S$sdneg
  }

  par(mfcol=c(4,1),oma=c(0.5,3,3,0))
#  par(oma=c(0,5,0,0))

  scorename<-c("log","quad","binary","class")
  algos<-c("BL-GL","BL-MoG","GH","LiNGAM","PC")
  order<-c(3,4,1,2)

  for ( scoreindex in 1:4 ) {
    score<-order[scoreindex]
    par(mfg=c(scoreindex,1))
    H2<-H[score,,]
    #H2[is.infinite(H2)]<-NA
    A1<-Ineg[score,,]
    A2<-Ipos[score,,]
    ylim<-c(min(min(A1[is.finite(A1)]),min(A2[is.finite(A2)]),min(H2[is.finite(H2)]))-0.00,
            max(max(A1[is.finite(A1)]),max(A2[is.finite(A2)]),max(H2[is.finite(H2)])))
    if (score == 1 ) {
      ylim[2]<-ylim[2]+0.3
    }

    B2<-array(0,c(nrow(A1),ncol(A1)))
    B1<-array(0,c(nrow(A1),ncol(A1)))

    B2[is.infinite(H2)]<-0.2
    B1[is.infinite(H2)]<-(-0.2)


    H2[is.infinite(H2)]<-ylim[2]-0.2
    A1[is.infinite(A1)]<-NA
    A2[is.infinite(A2)]<-NA

    print(ylim)
  
   par(mai=c(0,0.5,0,0))
    if ( score == 4 ) {
 
    } else {
    }#ylab=scorename[score],
    if ( score == 1 ) {
      B<-barplot(H2,beside=TRUE,ylim=c(ylim[1],3.2), space=c(0, 1),axes=FALSE,col=gray(seq(from=0.99,to=0.3,by=-0.7/5)))
      axis(side=2,outer=FALSE,labels=c("0","1.5","3"  ),at=c(0,1.5,3),cex.axis=1.8)
    } else if ( score == 2 ) {
      B<-barplot(H2,beside=TRUE,ylim=c(ylim[1],2.2 ), space=c(0, 1),axes=FALSE,col=gray(seq(from=0.99,to=0.3,by=-0.7/5)))
      axis(side=2,outer=FALSE,labels=c("0","1","2"  ),at=c(0,1,2),cex.axis=1.8)
    } else if ( score == 3 ) {
      B<-barplot(H2,beside=TRUE,ylim=c(ylim[1],1.1  ), space=c(0, 1),axes=FALSE,col=gray(seq(from=0.99,to=0.3,by=-0.7/5)))
      axis(side=2,outer=FALSE,labels=c("0","1/2","1"  ),at=c(0,0.5,1),cex.axis=1.8)
    } else if ( score == 4 ) {
       B<-barplot(H2,beside=TRUE,ylim=c(ylim[1],max(1.1,ylim[2])  ),
                axes=FALSE,col=gray(seq(from=0.99,to=0.3,by=-0.7/5)))#ylab=scorename[score],
 
      legend(x=55.0,y=1.0,c("BL-GL","BL-MoG","GH","LiNGAM","PC"),bty='n',fill=gray(seq(from=0.99,to=0.3,by=-0.7/5)),cex=1.5)
      axis(side=2,outer=FALSE,labels=c("0","1/2","1"  ),at=c(0,0.5,1),cex.axis=1.8)
    }
    #print(A1)
    #print(H2)
    #print(A2)
    arrows(B,A1,B,A2,code=3,angle=90,length=0.05)
    arrows(B,H2+B1,B,H2+B2,code=2,angle=45,length=0.15)
  }
   mtext("binary",side=2,at=7/8,outer=TRUE,cex=1.5,line=1)
   mtext("class",side=2,at=5/8,outer=TRUE,cex=1.5,line=1)
   mtext("log",side=2,at=3/8,outer=TRUE,cex=1.5,line=1)
   mtext("quadratic",side=2,at=1/8,outer=TRUE,cex=1.5,line=1)

   offset<-0.0525
   scale<-0.92
   mtext("Set 1",side=3,at=offset+scale*1/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 2",side=3,at=offset+scale*3/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 3",side=3,at=offset+scale*5/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 4",side=3,at=offset+scale*7/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 5",side=3,at=offset+scale*9/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 6",side=3,at=offset+scale*11/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 7",side=3,at=offset+scale*13/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 8",side=3,at=offset+scale*15/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 9",side=3,at=offset+scale*17/20,outer=TRUE,cex=1.5,line=1)
   mtext("Set 10",side=3,at=offset+scale*19/20,outer=TRUE,cex=1.5,line=1)



}
#use pdf(file='pseudorealplot.pdf',width=16.5,height=4/5*9/13*16.5)   its inches but a good picture when just give the centimeters
