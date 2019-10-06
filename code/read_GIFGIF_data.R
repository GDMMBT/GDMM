file=read.csv(file.choose(),header = T)
attach(file)

combine=c(as.character(left),as.character(right))
combine=as.factor(combine)
combine
w=matrix(0,length(levels(combine)),length(levels(combine)))

combine=as.numeric(combine)
left=combine[1:length(left)]
right=combine[(length(left)+1):length(combine)]

choice=as.numeric(choice)
file1=data.frame(left,right,choice)

###building winning matrix w_i,j denote the number of times that i beats j.

for(i in 1:dim(file1)[1]){
  if(file1[,3][i]==1){w[file1[i,1],file1[i,2]]=w[file1[i,1],file1[i,2]]+1}
  else if(file1[,3][i]==3){w[file1[i,2],file1[i,1]]=w[file1[i,2],file1[i,1]]+1}
}

n=dim(w)[1]
n

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

###maxinum number of paired comparison
m.max=max(diag(D))
m.max

###The laplacian matirx
L=D-A

###The second smallest eigenvalue
aaa=sort(eigen(L)$values)
aaa[2]
