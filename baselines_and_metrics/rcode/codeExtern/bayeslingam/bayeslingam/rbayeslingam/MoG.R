#Includes the functions for the calculation of the loglikelihood of the MoG-model
#first learn how to change from vector to matrix and vice versa

#H[lower.tri(H,diag=TRUE)] returns the vector form in the order:
#> H
#     [,1] [,2] [,3]
#[1,]    1    0    0
#[2,]    2    4    0
#[3,]    3    5    6


#WITHIN THESE WE NEED TO MAINTAIN THE "DATAGRADIENT" AND "DATAHESSIAN"
#SO THAT FINALLY WHEN x=y-X*b, df/db=(df/dx)*(dx/db)
#AND d2f/db2=(d2f/dxdx)*(dx/db)^2+(df/dx)*(d2x/db2)

vec2matrix<-function(h) {
  #the dimension of the matrix
  #d(d+1)/2 = length(h)
  d<-(sqrt(1+8*length(h))-1)/2

  H<-array(0,c(d,d))
  H[lower.tri(H,diag=TRUE)]<-h
  H<-H+t(H)
  diag(H)<-diag(H)/2  

  H
}

value<-function(val) {
  res<-val                #stupid function to get rid of attributes
  attributes(res)<-NULL   #in a situation like attr(b,"k")<-a*attr(c,"k")
                          #attributes of become attributes of attibute
  res
}

#REMEMBER TO TREAT THE FIRST FUNCTIONS AS JUST VECTORIZED FORM FOR A PARTICULAR SCALAR X-VALUE!!!

#for verbal
DPL<-0

logGaussian<-function(mu,sigmasq,x,derivatives=TRUE) {

  r<-x-mu
  rsq<-r^2

  val<-(-1/2)*log(2*pi)-1/2*log(sigmasq)-1/2*rsq/sigmasq

  if ( derivatives ) {
    #dL/dmu,dL/dsigmasq
    attr(val,"gradient")<-cbind(r/sigmasq,-1/(2*sigmasq)+1/2*rsq/(sigmasq^2))

    #dL/dx
    attr(val,"datagradient")<-(-r/sigmasq)

    #dL2/dmu2,dL2/dmudsigmasq,dL2/dsigmasq2
    attr(val,"hessian")<-cbind(-1/sigmasq, -r/(sigmasq^2),
                               1/(2*(sigmasq^2))-rsq/sigmasq^3)

    #dL2/dmudx,dL2/dsigmasqdx,dL2/dx2
    attr(val,"datahessian")<-cbind(1/sigmasq,r/(sigmasq^2),-1/sigmasq)
  }

  val
}

#all column products in order: 11,12,13,1N,22,23,2N,33,...

colProduct<-function(G1,G2=c()) {
  #we want only products with one column from each matrix
  G1<-as.matrix(G1)

  if ( is.null(G2) ) {
    #if G2 is not defined then G1 is used for G2
    N<-ncol(G1)
    R<-array(NA,c(nrow(G1),(N+1)*N/2))
    k<-1
    for (i in 1:N) {
      for (j in i:N) { #starting from i
        R[,k]<-G1[,i]*G1[,j]
        k<-k+1
      }
    }

  } else {
    G2<-as.matrix(G2)
    R<-array(NA,c(nrow(G1),ncol(G1)*ncol(G2)))
    k<-1
    for (i in 1:ncol(G1)) {
      for (j in 1:ncol(G2)) {
        R[,k]<-G1[,i]*G2[,j]
        k<-k+1
      }
    }
  }

  R
}

gaussian<-function(mu,sigmasq,x,derivatives=TRUE) {

  #exponentiating from logGaussian
  val<-logGaussian(mu,sigmasq,x,derivatives)

  res<-exp(value(val))

  if ( derivatives ) {
    #dG/dmu,dG/dsigmasq
    attr(res,"gradient")<-value(res)*attr(val,"gradient")
    #multiplication goes row by row

    #dG/dx
    attr(res,"datagradient")<-value(res)*attr(val,"datagradient")
    #also row by row multiplication

    #d(dG/dt)/ds=d(exp(L)*dL/dt)/ds=exp(L)*(dL/ds*dL/dt+d2L/dtds)
    #dL2/dmu2,dL2/dmudsigmasq,dL2/dsigmasq2
    #->dG2/dmu2,dG2/dmudsigmasq,dG2/dsigmasq2

    attr(res,"hessian")<-value(res)*(colProduct( attr(val,"gradient") )+
                                      attr(val,"hessian") )

    dsdt<-colProduct(cbind( attr(val,"gradient"),attr(val,"datagradient") ),
                           attr(val,"datagradient") )

    #dL2/dmudx,dL2/dsigmasqdx,dL2/dx2->dL2/dmudx,dL2/dsigmasqdx,dL2/dx2
    attr(res,"datahessian")<-value(res)*(dsdt+attr(val,"datahessian"))
  }

  res
}

piiGaussian<-function(pii,mu,sigmasq,x,derivatives=TRUE) {

  #this just multiplies with pii the value from gaussian
  val<-gaussian(mu,sigmasq,x,derivatives)

  res<-pii*value(val)

  if ( derivatives ) {

    #dP/dt=pii*dG/dt or dP/dpii=G
    attr(res,"gradient")<-cbind(val,pii*attr(val,"gradient"))

    attr(res,"datagradient")<-pii*attr(val,"datagradient")

    #with the hessian there is more trouble

    #dP2/dpii2=0,dP2/dpiidt=dG/dt, rest as before
    attr(res,"hessian")<-cbind(0,attr(val,"gradient"),pii*attr(val,"hessian"))

    #dP2/dpidx=dG/dx, rest as before
    attr(res,"datahessian")<-cbind(attr(val,"datagradient"),
                                   pii*attr(val,"datahessian"))
  }
  res
}

mixtureGaussian<-function(p,x, derivatives=TRUE) {

  N<-length(x) #the length of data
  M<-length(p)/3 #the number of mixtures

  # the variables decoded from the parameter vector
  pii<-p[seq(from=1,to=(3*M),by=3)]
  mu<-p[seq(from=2,to=(3*M),by=3)]
  sigmasq<-p[seq(from=3,to=(3*M),by=3)]

  res<-rep(0,N)
  G<-array(0,c(N,3*M)) #each mixture has 3 parmeters
  DG<-rep(0,N)

  m<-3*M #total number of parameters
  H<-array(0,c(N,m*(m+1)/2)) #this many parameter-parameter elements

  #in addition the datahessian elements,
  #in total m+1 columns!
  DH<-array(0,c(N,m+1))


  s<-1 #for hessian indexing
  for(g in (1:M)) {
    val<-piiGaussian(pii[g],mu[g],sigmasq[g],x,derivatives)

    #value is just added together
    res<-res+value(val)

    if ( derivatives ) {
      #insert the three columns of the gradient to the right place
      #because dM/dp_i=d/dp_i sum_j(P(p_j,x))=dP(p_i,x)/dp_i
      index<-((g-1)*3+1):(g*3)#for gradient and datahessian indexing
      G[,index]<-attr(val,"gradient")

      #dM/dx=sum_i(dP(p_i),x)/dx)
      DG<-DG+attr(val,"datagradient")

      #first of all d2M/dp_idp_j=0! (parameters of two different gaussians)
      #and d2M/dp_idp_j=d2P(p_i,p_j,x)/dp_idp_j /
      # this is the case when to parameters are in the same gaussian
      #H=[dP(g=1),3*(M-g) rows of zeros,dP(g=2),3*(M-g) rows of zeros,...]
      #

      #fill in the for parameters within same gaussian
      H[,s:(s+2)]<-attr(val,"hessian")[,1:3]
      s<-s+2+1


      s<-s+3*(M-g)

      H[,s:(s+1)]<-attr(val,"hessian")[,4:5]
      s<-s+1+1

      s<-s+3*(M-g)

      H[,s]<-attr(val,"hessian")[,6]
      s<-s+1

      s<-s+3*(M-g)

      #d2M/dp_jdx=d/dx(dP(p_j,x)/dp_j)=d2P(p_j,x)/dp_jdx
      #copying to the right place
      DH[,index]<-attr(val,"datahessian")[,1:3]

      #d2M/dx2=sum_i d2P/dx2
      DH[,m+1]<-DH[,m+1]+attr(val,"datahessian")[,4]
      #the fourth component is d2P/dx2
    }# if derivatives
  }#for

  if ( derivatives ) {
    attr(res,"gradient")<-G
    attr(res,"hessian")<-H
    attr(res,"datagradient")<-DG
    attr(res,"datahessian")<-DH
  }

  res
}

logMixtureGaussian<-function(p,x,derivatives=TRUE) {

  val<-mixtureGaussian(p,x,derivatives)

  res<-log(value(val))
  #L=log(M)
  if ( derivatives ) {
    #dL/dt=1/M*dM/dt

    attr(res,"gradient")<-attr(val,"gradient")/value(val)#division row by row
    attr(res,"datagradient")<-attr(val,"datagradient")/value(val)

    #d2L/dtds=d/ds(1/M*dM/dt)=d/ds(1/M)*dM/dt+1/M*(d2M/dtds)
    #=-1/M^2*(dM/ds)*(dM/dt)+1/M*(d2M/dtds)
    #notice that these produces nonzeros elements 
    #for parameters from different gaussians!

    attr(res,"hessian")<-(-colProduct(attr(val,"gradient"))/(value(val)^2)+
                          attr(val,"hessian")/value(val)) 

    dsdt<-colProduct(cbind(attr(val,"gradient"),attr(val,"datagradient")),
                     attr(val,"datagradient"))

    attr(res,"datahessian")<-(-dsdt/(value(val)^2)+
                             attr(val,"datahessian")/value(val))
  }
  res
}

#HERE WE ARE NOT DEALING WITH JUST A VECTORIZED FORM SINCE THE VALUES
#FOR DIFFERENT DATA ARE ADDED TOGETHER!

sumLogMixtureGaussian<-function(p,x,derivatives=TRUE) {
  val<-logMixtureGaussian(p,x,derivatives)

  res<-sum(value(val))
  if ( derivatives ) {
    #leaving dimension 2 = columns as it is
    #dS/dp=sum_i dP(x_i,p)/dp

    attr(res,"gradient")<-c(apply(attr(val,"gradient"),2,sum)) 
    #now the gradient is a vector!

    #for data gradient just copy
    #dS/dx_j=(d/dx_j) sum_i P(x_i,p) = dP(x_j,p)/dx_j
    attr(res,"datagradient")<-attr(val,"datagradient")

    #d2S/dp1dp2=(d/dp2) (sum_i dP(x_i,p1,p2)/dp1)=sum_i dP(x_i,p1,p2)/dp1dp2

    h<-apply(attr(val,"hessian"),2,sum) #change the hessian to its matrix form!

    attr(res,"hessian")<-vec2matrix(h)

    #d2S/dpdx_j=(d/dp) dP(x_j,p)/dx_j=d2P(x_j,p)/dx_jdp
    #d2S/dx_jdx_j=(d/dx_j) dP(x_j,p)/dx_j=d2P(x_j,p)/dx_jdx_j
    #cross terms from data would be zero, so not in here!
    attr(res,"datahessian")<-attr(val,"datahessian")
  }

  res
}

parchange<-function(p,print.level=DPL,derivatives=TRUE) {
  #3M parameters to 3M parameters for now!
  m<-length(p)
  M<-m/3

  alphaindex<-seq(from=1,to=m,by=3)
  muindex<-seq(from=2,to=m,by=3)
  gammaindex<-seq(from=3,to=m,by=3)

  alpha<-p[alphaindex]
  mu<-p[muindex]
  gamma<-p[gammaindex]

  q<-rep(NA,m)

  ea<-exp(alpha)
  sa<-sum(ea)
  pii<-ea/sa

  sigmasq<-exp(-gamma)

  q[alphaindex]<-pii #parameter change
  q[muindex]<-mu
  q[gammaindex]<-sigmasq

  if ( derivatives ) {
    #print(q)
    #because this a vectorized function the gradient is a 
    #m functions on m parameters
    G<-array(0,c(m,m))

    #df_i/dmu_j=1 if f_i is a mu transform and i=j, otherwise zero 
    G[cbind(muindex,muindex)]<-1

    #df_i/dgamma_j=-exp(-gamma) if f_i is a gamma transform and i=j,
    # otherwise zero 
    G[cbind(gammaindex,gammaindex)]<-(-1)*sigmasq

    #dpi_i/dalpha_i=pi_i-pi_i^2
    #dpi_i/dalpha_k=-pi_i*pi_k
    #grad pi = diag(pi)-pi%*%t(pi)
    G[alphaindex,alphaindex]<-diag(pii)-pii%*%t(pii)

    #print(G)
    #because the functions were linear all second derivatives are zero!
    H<-array(0,c(m,m,m))#Hs are
    #let the first index be the fi index! 

    #H[muindex,,]<-0!
    #print(H)
    #d2F/dgamma2=exp(-gamma)
    #so for gamma transforms hessian is has 
    #exp(-gamma) on one element on the diagonal

    H[cbind(gammaindex,gammaindex,gammaindex)]<-sigmasq

    for ( i in 1:M ) {
      k<-alphaindex[i]
      H[k,alphaindex,alphaindex]<-2*pii[i]*pii%*%t(pii)

      H[k,k,alphaindex]<-H[k,k,alphaindex]-pii[i]*t(pii)

      H[k,alphaindex,k]<-H[k,alphaindex,k]-pii[i]*t(pii)

      if ( M > 1 ) {
        diag(H[k,alphaindex,alphaindex])<-diag(H[k,alphaindex,alphaindex])+
                                          (-pii[i])*pii
      } else {
        H[k,alphaindex,alphaindex]<-H[k,alphaindex,alphaindex]+(-pii[i])*pii
      }
      H[k,k,k]<-H[k,k,k]+pii[i]
    }

    attr(q,"gradient")<-G
    attr(q,"hessian")<-H
  }#if derivatives

  q
}

parSumLogMixtureGaussian<-function(p,x,derivatives=TRUE) {
  #this should work with the changed parameters
  q<-parchange(p,derivatives)

  #Hq is a 3-D matrix
  #attributes(q)<-NULL

  val<-sumLogMixtureGaussian(q,x,derivatives)

  res<-value(val)

  if ( derivatives ) {
    #G[i]=sum_j Gf[j]*Gy[j,i] = sum_j t(Gy)[i,j]*Gf[j]=(t(Gy)%*%Gf)[i]
    #G=t(Gy)%*%Gf
    G<-t(attr(q,"gradient"))%*%attr(val,"gradient")

    Hq<-attr(q,"hessian")
    Gf<-attr(val,"gradient")

    H2<-array(0,c(length(q),length(q)))
    for (j in 1:(dim(H2)[1])) {
      H2<-H2+Gf[j]*Hq[j,,]
    }

    H<-t(attr(q,"gradient"))%*%attr(val,"hessian")%*%attr(q,"gradient")+H2

    #res<-(-1)*res
    attr(res,"gradient")<-c(G)
    attr(res,"hessian")<-H

    #now implement the change in the data gradient
    #first of all d2f/dxdgamma=d2f/dxdsigma*dsigma/dgamma 
    DG<-attr(val,"datagradient")
    DH<-attr(val,"datahessian")

    DH[,1:(ncol(DH)-1)]<-DH[,1:(ncol(DH)-1)]%*%t(attr(q,"gradient"))
    attr(res,"datagradient")<-attr(val,"datagradient")
    #here the transform is not yet applied!
    attr(res,"datahessian")<-DH
  }

 res
}



residual<-function(b,y,X,derivatives=TRUE) {
  #X is always a matrix

  #hessian not needed
  #anyway full of zeros

  if ( ncol(X)==1) {
    r<-y-b*X
    if ( derivatives ) {
      attr(r,"gradient")<-(-X)
      #attr(r,"hessian")<-array(0,c(1,1,1))
    }
  } else if ( ncol(X) == 0 ) {
    r<-y
    if ( derivatives ) {
      attr(r,"gradient")<-array(0,c(nrow(X),0))
      #attr(r,"hessian")<-array(0,c(0,0,0))#????
    }
  } else {
    r<-y-X%*%b
    if ( derivatives ) {
      attr(r,"gradient")<-(-X)
      #attr(r,"hessian")<-array(0,c(length(b),length(b),length(b)))
    }
  }

  r
}

 bParSumLogMixtureGaussian<-function(p,y,X,derivatives=TRUE) {
  #this function return the likelihood needed 
  #for linear mixture gaussian network modelling

  bs<-ncol(X)#X should always be a matrix!!!!

  ps<-length(p) #parameters

  b<-p[(ps-bs+1):ps]

  r<-residual(b,y,X,derivatives)
  if ( derivatives ) {
    rG<-attr(r,"gradient")
  }

  #The formulas below apply only in the case where r is linear for b!!!!!!!!!!!
  #so hessian(r) is zeros
  #
  #rH<-attr(r,"hessian")

  val<-parSumLogMixtureGaussian(p[1:(ps-bs)],value(r),derivatives)
  res<-value(val)

  if ( derivatives ) {
   mG<-attr(val,"gradient")
   mH<-attr(val,"hessian")
   mdG<-attr(val,"datagradient")
   mdH<-attr(val,"datahessian")

   G<-c(mG,mdG%*%rG)

   H<-array(NA,c( length(p),length(p) ) )

   H[1:(ps-bs),1:(ps-bs)]<-mH

    if ( bs != 0 ) { #only if there are bs!!!!
      if (length(rG) != 1) {
        H[1:(ps-bs),(ps-bs+1):ps]<-t(as.matrix(mdH[,1:(ncol(mdH)-1)]))%*%
                                   as.matrix(rG)
        H[(ps-bs+1):ps,1:(ps-bs)]<-t(H[1:(ps-bs),(ps-bs+1):ps])
      } else {
        H[(ps-bs+1):ps,1:(ps-bs)]<-mdH[,1:(ncol(mdH)-1)]*rG
        H[1:(ps-bs),(ps-bs+1):ps]<-mdH[,1:(ncol(mdH)-1)]*rG
      }

      hb<-c(apply(colProduct(as.matrix(rG))*mdH[,ncol(mdH)],2,sum))

      H[(ps-bs+1):ps,(ps-bs+1):ps]<-vec2matrix(hb)
    }

    attr(res,"gradient")<-G
    attr(res,"hessian")<-H
  }#if derivatives

  res
}













