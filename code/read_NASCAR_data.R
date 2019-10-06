NASCAR=read.table(file.choose(),header = T)
attach(NASCAR)
n=83
nrace=36
init=rep(0,n)
tiny=10^(-4)

#k is a list of total number of drivers for each race, there are 36 races.
#k[2:37] records total number of drivers for race 1 to race 36.
k=rep(0,(nrace+1))
for(i in 1:36){
  k[i+1]=sum(NASCAR$Race==i)}

#playersmatrix with each row represents the players ID with ranking orders.
playersmatrix=matrix(0,36,43)
for(i in 1:36){
  playersmatrix[i,(1:k[i+1])]=DriverID[(sum(k[1:i])+1):(sum(k[1:i])+k[i+1])]
}


###A_ij is the total number of comparisons that item i and item j are compared.
A=matrix(0,83,83)
for(k in 1:36){
  for(i in 1:83){
    if(sum(playersmatrix[k,]==i)==1){A[i,playersmatrix[k,]]=A[i,playersmatrix[k,]]+1}
  }
}

###D_ii is the total number of paried comparisons that i item has.
D=matrix(0,n,n)
for(i in 1:n){
  D[i,i]=sum(A[i,])
}

###maxinum number of paired comparison
m.max=max(diag(D))
m.max

###The laplacian matirx
L=D-A

###The second smallest eigenvalue
aaa=sort(eigen(L)$values)
aaa[2]


