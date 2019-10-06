chess=read.csv(file.choose(),header = T)
attach(chess)
n=max(White.Player..,Black.Player..)

###Winning matrix. w_ij is the number of observed paired comparison of items i and j such that i wins j.
w=matrix(0,n,n)
for(i in 1:dim(chess)[1]){
  if(chess[,4][i]==1){w[chess[i,2],chess[i,3]]=w[chess[i,2],chess[i,3]]+1}
  else if(chess[,4][i]==0){w[chess[i,3],chess[i,2]]=w[chess[i,3],chess[i,2]]+1}
}

###tie matrix. t_ij is the number of observed paired comparison of items i and j such that there is a tie.
t=matrix(0,n,n)
for(i in 1:dim(chess)[1]){
  if(chess[,4][i]==0.5){t[chess[i,2],chess[i,3]]=t[chess[i,2],chess[i,3]]+1}
}
t=t+t(t)

###comparison matrix with s_ij as the number of observed paired comparisons of items i and j such that either i wins j or there is tie outcome.
s=w+t

###A_ij is the total number of comparisons that item i and item j are compared.
A=matrix(0,n,n)
for(i in 1:n){
  for(j in 1:n){
    if((w[i,j]+w[j,i])!= 0){A[i,j]=w[i,j]+w[j,i]}
  }
}

###D_ii is the total number of paried comparisons that i item has.
D=matrix(0,n,n)
for(i in 1:n){
  D[i,i]=sum(A[i,])
}

m.max=max(D)

###The laplacian matirx
L=D-A

###The second smallest eigenvalue
aaa=sort(eigen(L)$values)
aaa[2]
