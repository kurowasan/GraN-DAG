sampleSize <- 500
x <- matrix(0,sampleSize,3)

for(i in 1:sampleSize)
{
  a <- sample(4,1)
  x[i,1:3] <-  switch(as.character(a),
         "1"= c(0,0,0),
         "2"= c(0,1,1),
         "3"= c(1,0,1),
         "4"= c(1,1,0))
}
show(mutual_hsic(x[,1],x[,2],x[,3]))

show(indtest_hsic(x[,2],x[,3]))
show(indtest_hsic(x[,1],x[,3]))
show(indtest_hsic(x[,2],x[,1]))
