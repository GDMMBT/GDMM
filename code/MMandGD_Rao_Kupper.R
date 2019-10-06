###MM algorithm for Rao_Kupper. c=theta
###s is comparison matrix with s_ij as the number of observed paired comparisons of items i and j such that either i wins j or there is tie outcome.
Bayesian.MM=function(init,tiny,alpha,beta,s,c){
  n=dim(s)[1]
  iterationmatrix=matrix(0,100000,n)
  t=0
  x1=init
  x0=x1+rep(10*tiny,n)
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(i in 1:n){
      u=c()
      l=c()
      for(j in 1:n){
        if(j != i){
          u[j]=s[i,j]
          l[j]=((s[i,j])/(exp(x0[i])+c*exp(x0[j])))+((c*s[j,i])/(exp(x0[j])+c*exp(x0[i])))
        }
        else{u[j]=0
        l[j]=0
        }
      }
      x1[i]=log((alpha-1+sum(u))/(beta+sum(l)))
    }
    t=t+1
    iterationmatrix[t,]=x1
  }
  iterationmatrix=iterationmatrix[1:t,1:n]
  return(list(Iterationmatrix=iterationmatrix))
}

#####################################################################################################################
#value settings for experiments MM:
#only for sample case
#s=s0
#n=n0

init=rep(0,n)
tiny=10^(-4)
c=sqrt(2)

ImatrixMM=Bayesian.MM(init,tiny,11,10,s,c)$Iterationmatrix
ImatrixMM.21=Bayesian.MM(init,tiny,2,1,s,c)$Iterationmatrix
ImatrixMM.1=Bayesian.MM(init,tiny,1.1,0.1,s,c)$Iterationmatrix
ImatrixMM.01=Bayesian.MM(init,tiny,1.01,0.01,s,c)$Iterationmatrix
ImatrixMM0=Bayesian.MM(init,tiny,1,0,s,c)$Iterationmatrix

dim(ImatrixMM)
dim(ImatrixMM.21)
dim(ImatrixMM.1)
dim(ImatrixMM.01)

dim(ImatrixMM0)

######################################################################################################################
###GD algorithm for Rao_Kupper. c=theta
###s is comparison matrix with s_ij as the number of observed paired comparisons of items i and j such that either i wins j or there is tie outcome.
Bayesian.GD=function(init,tiny,alpha,beta,s,m.max,c){
  d=c()
  dx=c()
  n=dim(s)[1]
  iterationmatrix=matrix(0,100000,n)
  t=0
  x1=init
  x0=x1+rep(10*tiny,n)
  while(sum(abs(x0-x1)>tiny)>0) {
    x0=x1
    for(i in 1:n){
      for(j in 1:n){
        if (j != i) {d[j]=(s[i,j]*(1-1/(1+c*exp(x0[j]-x0[i]))))-s[j,i]/(1+exp(x0[j]-x0[i])/c)}
        else{d[j]=0}
      }
      dx[i]=sum(d)+(alpha-1)-beta*exp(x0[i])
    }             #this returns dx
    x1=x0+(1/((m.max/2)+beta))*dx
    t=t+1
    iterationmatrix[t,]=x1
  }
  iterationmatrix=iterationmatrix[1:t,1:n] #total rows minus one = number of iterations
  return(list(Iterationmatrix=iterationmatrix))
}


#value settings for experiments GD:
#only for sample case
#s=s0
#n=n0 
#m.max=mm.max

init=rep(0,n)
tiny=10^(-4)
c=sqrt(2)

ImatrixGD=Bayesian.GD(init,tiny,11,10,s,m.max,c)$Iterationmatrix
ImatrixGD.21=Bayesian.GD(init,tiny,2,1,s,m.max,c)$Iterationmatrix
ImatrixGD.1=Bayesian.GD(init,tiny,1.1,0.1,s,m.max,c)$Iterationmatrix
ImatrixGD.01=Bayesian.GD(init,tiny,1.01,0.01,s,m.max,c)$Iterationmatrix

ImatrixGD0=Bayesian.GD(init,tiny,1,0,s,m.max,c)$Iterationmatrix

dim(ImatrixGD)
dim(ImatrixGD.21)
dim(ImatrixGD.1)
dim(ImatrixGD.01)

dim(ImatrixGD0)