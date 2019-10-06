Bayesian.GD=function(init,tiny,alpha,beta,w,m.max){
  d=c()
  dx=c()
  n=dim(w)[1]
  iterationmatrix=matrix(0,100000,n)
  i=0
  x1=init
  x0=x1+rep(10*tiny,length(x1))
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(j in 1:n){
      for(k in 1:n){
        if (k != j) {d[k]=w[j,k]-(w[j,k]+w[k,j])/(1+exp(x0[k]-x0[j]))}
        else{d[k]=0}
      }
      dx[j]=sum(d)+(alpha-1)-beta*exp(x0[j])
    }             #this returns dx
    x1=x0+(1/((m.max/2)+beta))*dx
    i=i+1
    iterationmatrix[i,]=x1
  }
  iterationmatrix=iterationmatrix[1:i,1:n] #total rows minus one = number of iterations
  return(list(Iterationmatrix=iterationmatrix))
}

##################################################################################
Bayesian.MM=function(init,tiny,alpha,beta,w){
  n=dim(w)[1]
  iterationmatrix=matrix(0,100000,n)
  i=0
  x1=init
  x0=x1+rep(10*tiny,length(x1))
  u=c()
  l=c()
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(j in 1:n){
      for(k in 1:n){
        if(k != j){
          u[k]=w[j,k]
          l[k]=(w[j,k]+w[k,j])/(exp(x0[j])+exp(x0[k]))
        }
        else{u[k]=0
        l[k]=0
        }
      }
      x1[j]=log((alpha-1+sum(u))/(beta+sum(l)))
    }
    i=i+1
    iterationmatrix[i,]=x1
  }
  iterationmatrix=iterationmatrix[1:i,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

##################################################################################
init=rep(0,n)
tiny=10^-4

#value settings for experiments GD:
#only for sample case
#s=s0
#n=n0 
#m.max=mm.max

ImatrixGD=Bayesian.GD(init,tiny,11,10,w,m.max)$Iterationmatrix

ImatrixGD.21=Bayesian.GD(init,tiny,2,1,w,m.max)$Iterationmatrix
ImatrixGD.1=Bayesian.GD(init,tiny,1.1,0.1,w,m.max)$Iterationmatrix
ImatrixGD.01=Bayesian.GD(init,tiny,1.01,0.01,w,m.max)$Iterationmatrix
ImatrixGD0=Bayesian.GD(init,tiny,1,0,w,m.max)$Iterationmatrix


dim(ImatrixGD)
dim(ImatrixGD.21)
dim(ImatrixGD.1)
dim(ImatrixGD.01)
dim(ImatrixGD0)

##################################################################################

ImatrixMM=Bayesian.MM(init,tiny,11,10,w)$Iterationmatrix
ImatrixMM.21=Bayesian.MM(init,tiny,2,1,w)$Iterationmatrix
ImatrixMM.1=Bayesian.MM(init,tiny,1.1,0.1,w)$Iterationmatrix
ImatrixMM.01=Bayesian.MM(init,tiny,1.01,0.01,w)$Iterationmatrix
ImatrixMM0=Bayesian.MM(init,tiny,1,0,w)$Iterationmatrix

dim(ImatrixMM)
dim(ImatrixMM.21)
dim(ImatrixMM.1)
dim(ImatrixMM.01)
dim(ImatrixMM0)