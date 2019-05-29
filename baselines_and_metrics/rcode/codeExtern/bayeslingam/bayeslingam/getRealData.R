realdatadir<-sprintf('%s%s',cauzality_path,'/trunk/bayeslingam/data/realdata/')

getRealData<-function( i,verbal=2 ) {
  #This function returns real data from data set i
  X<-NULL
  if ( i == 1  ) { #first 8 data sets are from the causal challenge
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 1 data: Altitude, Temperature\n')
      cat('Proposed DAG: A->T\n')
    }
  } else if ( i== 2 ) {
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 2 data: Altitude, Precipitation\n')
      cat('Proposed DAG: A->P\n')
    }
  } else if ( i== 3 ) {
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 3 data: Longitude, Temperature\n')
      cat('Proposed DAG: L->T\n')
    }
  } else if ( i== 4 ) {
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 4 data: Sunshine, Altitude\n')
      cat('Proposed DAG: S<-A\n')
    }
  } else if ( i== 5 ) {
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 5 data: Abalone Length, Age \n')
      cat('Proposed DAG: L<-A\n')
    }
  } else if ( i== 6 ) {
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 6 data: Abalone Age, Weight\n')
      cat('Proposed DAG: A->W\n')
    }
  } else if ( i== 7 ) {
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 7 data: Abalone Diameter, Age\n')
      cat('Proposed DAG: D<-A\n')
    }
  } else if ( i== 8 ) {
    file<-sprintf('%spairs0%i.txt',realdatadir,i)
    X<-read.table(file,sep=" ")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge 8 data: Age, (hourly) Wages\n')
      cat('Proposed DAG: A->W\n')
    }
  } else if ( i== 9 ) {
    X<-read.table(sprintf('%s/abalone.data',realdatadir),sep=",")
    X<-X[,c(2,5,9)] #2nd is the length, 9th is the length
    names(X)<-c("length","weight","age")
    if ( verbal >= 2 ) {
      cat('Loaded Challenge extended data: Abalone Length, Weight, Age.\n')
      cat('Proposed DAG: A->L, (L->W), A->W\n')
    }
  } else if ( i== 10 ) {
    X<-read.table(sprintf('%s/economicindex.txt',realdatadir),sep=" ")
    names(X)<-c("Dow Jones","AAL")
    if ( verbal >= 2 ) {
      cat('Loaded economic index data:\nDow Jones index and Australian all ordinaries index at closing on 251 trading days ending 26 August 1994.\n')
      cat('Proposed DAG: D-A\n')
    }
  } else if ( i== 11 ) {
    X<-read.table(sprintf('%s/forex.txt',realdatadir),header=TRUE)
    X<-X[1:1000,c(2,3,8,9)]
    names(X)<-c("Aus","Brit","Jap","Swiss")
    if ( verbal >= 2 ) {
      cat('Loaded forex exchange rate data:\n Aus/US,Brit/US,Jap/US,SWIS/US\n')
      cat('Proposed DAG: (Aus-Brit-Jap),Aus->Swiss,Jap->Swiss \n')
    }
  } else if ( i== 12 ) {
    X<-read.table(sprintf('%s/nhanes3.dat',realdatadir),header=TRUE,na.strings=".")
    X<-X[1:2000,c(8,9,10,11,15)] #colestrol out 15, "Col"
    names(X)<-c("Weight","Height","Sys","Dia","Col")
    X<-X[!is.na(X[,1]) & !is.na(X[,2]) & !is.na(X[,3])& !is.na(X[,4]) &!is.na(X[,5]),]

     if ( verbal >= 2 ) {
       cat('Loaded health data:\n Weight,Height,Cholesterol,Systolic,Diastolic Blood Pressures\n')
       cat('Proposed DAG: Complicated \n')
     }
  } else if ( i== 13 ) {
    X<-getRealData(10)
    Y<-getRealData(12)
    n<-min(c(nrow(X),nrow(Y)))
    X<-data.frame(cbind(X[1:n,1],Y[1:n,4]))
    names(X)<-c("Dow Jones","Blood Pressure")
    if ( verbal >= 2 ) {
      cat('Loaded miex data:\n Dow Jones, Blood pressure\n')
      cat('Proposed DAG: Definitely Independent \n')
    }
  } else if ( i== 14 ) {
    X<-getRealData(11)
    set.seed(666)
    I<-sample(nrow(X),nrow(X))
    X[,1]<-X[I,1]
    X[,4]<-X[I,4]
    names(X)<-c("Aus","Brit","Jap","Swiss")
    if ( verbal >= 2 ) {
      cat('Loaded mixed forex exchange rate data:\n Aus/US,Brit/US,Jap/US,SWIS/US\n')
      cat('Proposed DAG: (Aus-Brit-Jap),Aus->Swiss,Jap->Swiss \n')
    }
  } else if ( i== 15 ) {
     X<-getRealData(9)
     set.seed(667)
     n<-nrow(X)
     X[,1]<-X[sample(n,n),1]
     X[,2]<-X[sample(n,n),2]
     X[,3]<-X[sample(n,n),3]
    names(X)<-c("length","weight","age")
    if ( verbal >= 2 ) {
      cat('Loaded shuffled data: Abalone Length, Weight, Age.\n')
      cat('Proposed DAG: Independent\n')
    }
  } else if ( i== 16 ) {
    X<-getRealData(12)
    Y<-getRealData(1)
    n<-min(c(nrow(X),nrow(Y)))
    X<-data.frame(cbind(X[1:n,c(1,2)],Y[1:n,1]))
    names(X)<-c("Height","Weight","Altitude")
     if ( verbal >= 2 ) {
       cat('Loaded mixed health data:\n Weight,Height,Altitude\n')
       cat('Proposed DAG: W-H , A\n')
     }
  } else if ( i== 17 ) {
    X<-read.table(sprintf('%s/arrhythmia.data',realdatadir),header=TRUE,na.strings=c("-","0.0"),sep=",")
    X<-X[,c(5,6,9)]
    X<-data.frame(X[apply(X<0.001,1,sum) == 0,])
    names(X)<-c("QRS","P-R" ,"P")
     if ( verbal >= 2 ) {
       cat('Loaded Arrhytmia data:\n QRS P-R P-interval\n')
       cat('Proposed DAG: QRS, PR-P\n')
     }

  } else if ( i== 18 ) {
    X<-read.table(sprintf('%s/forestfires.csv',realdatadir),sep=",",header=TRUE)
    X<-data.frame(X[,9:11])
    names(X)<-c("temp","RH","wind")
     if ( verbal >= 2 ) {
       cat('Loaded Forest fire data:\n TEMP, RH, WIND\n')
       cat('Proposed DAG: somehow connected\n')
     }

  } else if ( i== 19 ) {
    X<-getRealData(20)
    Y<-getRealData(18)
    n<-min(c(nrow(X),nrow(Y)))
    X<-data.frame(cbind(X[1:n,1],Y[1:n,2]))
     if ( verbal >= 2 ) {
       cat('Loaded Mixed data:\n Eruption time, Humidity\n')
       cat('Proposed DAG: independent\n')
     }
  } else if ( i== 20 ) {
    data(faithful)
    X<-faithful
    names(X)<-c("eruption","waiting")
    if ( verbal >= 2 ) {
      cat('Loaded Old Faithful Geysir data: Eruption, Waiting \n')
      cat('Proposed DAG: E->W\n')
    }
  } else {
    cat( 'Real data index not within range 1:20\n' );
  }
  X
}