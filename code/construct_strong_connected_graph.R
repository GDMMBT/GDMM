### Note that ###
### we use this code to extract the strongly connected sample graph from n sample nodes.
### we set n=1000 for chess
### n=285 for happiness (H)
### n=332 for contempt (C)
### n=350 for amusement (A)
### n=1000 for chess
### At the end, we will get strong connected sample graph for each case
### n0 is the number of nodes contained in each sample graph
### following the code below, we will get n0=251 for (A)
### n0=256 for (C)
### n0=252 for (A)
### n0=985 for chess
### w0 is the new winning matrix for the sample graph

library(igraph)
n=1000
cc=order(diag(D),decreasing = T)[1:n]
dd=diag(D)[cc]
DD=matrix(0,n,n)
diag(DD)=dd
cc1=sort(diag(D),decreasing = T)[1:n]

ww=matrix(0,n,n)
for(i in 1:n){
  for(j in 1:n){
    ww[i,j]=w[cc[i],cc[j]]}
}

graph1=graph_from_adjacency_matrix(ww,mode = "directed")
component=clusters(graph1,mode = "strong")
component

cc0=cc[which(component$membership == which.max(component$csize))]
n0=length(cc0)
n0

w0=matrix(0,n0,n0)
for(i in 1:n0){
  for(j in 1:n0){
    w0[i,j]=w[cc0[i],cc0[j]]}
}

#ONLY FOR CHESS DATA WITH TIES
s0=matrix(0,n0,n0)
for(i in 1:n0){
  for(j in 1:n0){
    s0[i,j]=s[cc0[i],cc0[j]]}
}

graph.1=graph_from_adjacency_matrix(w0,mode = "directed")
is.connected(graph.1,mode = "weak")
is.connected(graph.1,mode = "strong")

###AA_ij is the total number of comparisons that item i and item j are compared.
AA=matrix(0,n0,n0)
for(i in 1:n0){
  for(j in 1:n0){
    AA[i,j]=A[cc0[i],cc0[j]]}
}

###DD_ii is the total number of paried comparisons that i item has.
DD=matrix(0,n0,n0)
for(i in 1:n0){
  DD[i,i]=sum(AA[i,])
}

###Laplacian matrix
LL=DD-AA

###second smallest eigenvalue
sort(eigen(LL)$values)[2]

###maximum number of paired comparison
mm.max=max(diag(DD))
mm.max

###total number of comparisons
sum(w0)
