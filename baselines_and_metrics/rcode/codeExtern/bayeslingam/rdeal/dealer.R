#dealer, stub for the deal function

dealer<-function(X) {
  #input: a matrix with the data of each node as a column
    library('deal') #make sure library is in use
    #this follows the document: The deal Package February 16, 2008
    nvars<-ncol(X)
    names(X)<-paste("X",1:ncol(X),sep="") #names for dealer

    Xf<-data.frame(X) #make a data frame of X
    nw<-network(Xf)

    prior<-jointprior(nw,100,phiprior="heckerman")

    r <- list()

    if ( nvars <= 4 ) {
      #prior must be created with a count, otherwise errors in networkfamily!

      nf<-networkfamily(Xf,nw=nw,prior=prior)
      cat('done\n');
      #    print(nf);

      R<-nf$nw

      n<-length(R)

      tdags<-rep('',n)
      dim(tdags)<-c(n,1)
      r$loglike<-rep(NA,n)


      for (j in 1:n ) { #must gather the data into proper format
         tdags[j]<-modelstring(R[[j]])
         r$loglike[j]<-R[[j]]$score
      }



    } else {
      #use here the autosearch function somehow
       nw<-getnetwork(learn(nw=nw,df=Xf,prior=prior))
      cat('learn done\n');
       nf<-autosearch(initnw=nw,data=Xf,prior=prior,trace=FALSE,removecycles=TRUE)
      cat('done\n');

      tdags<-modelstring(nf$nw)
      dim(tdags)<-c(1,1)
      r$loglike<-nf$nw$score
    }

    r$prob <- exp(r$loglike - max(r$loglike))
    r$prob<-r$prob/sum(r$prob)
    r$DAGs<-bdag.to.cdag(tdag.to.bdag(tdags))
    #print(r)
  r
}

#functions for converting and comparing dags

equalDAG<-function( dealstring, dagvector ) {
  nvars<-max(dagvector)
  B1<-dagvector2B(dagvector)#DAG can be == only by the B-notation
  B2<-dealstring2B(dealstring,nvars)
  #print(B1)
  #print(B2)
  all(B1 == B2)
}

dealstring2B<-function(dealstring,nvars) {
  varnames<-paste("X",c(1:nvars),sep="")#the names of the variables in order!
  B<-array(0,c(nvars,nvars))
  s<-substr(dealstring,start=2,stop=nchar(dealstring)-1) 
  #first rip off the beginning [ and ending ]
  s<-strsplit(s,"\\]\\[")[[1]]#the split with ][

  for (i in 1:length(s)) {
    si<-s[i]
    k<-strsplit(si,"\\|")[[1]]
    node<-which(varnames==k[1])
    if ( length(k) > 1 ) {
      parentsstring<-k[2]
    } else {
      parentsstring<-""
    }
    parents<-c()
    while( nchar(parentsstring) != 0 ) {
      for (j in 1:length(varnames)) {
        if ( substr(parentsstring,start=1,
                    stop=nchar(varnames[j])) == varnames[j] ) {
          parents<-c(parents,j)
          parentsstring<-substr(parentsstring,start=nchar(varnames[j])+2,
                                stop=nchar(parentsstring))
        }
      }
    }
    B[node,parents]<-1
  }
  B
}
