library(Matrix)
linearsolve1=function(y,X){
  m=X
  m=crossprod(X,m)
  b=crossprod(X,y)
  return=solve(m)%*%b
}

linearsolve2=function(y,X){
  m=X
  m=crossprod(X,m)
  b=crossprod(X,y)
  lufact=expand(lu(m))
  l=as.matrix(lufact$L)
  u=as.matrix(lufact$U)
  p=lufact$P
  b=p%*%b
  for(i in 1:nrow(l)){
    db=b[i]/l[i,i]*as.vector(l[,i])
    db[i]=0
    b=b-db
    b[i]=b[i]/l[i,i]
  }
  for(i in rev(1:nrow(u))){
    db=b[i]/u[i,i]*as.vector(u[,i])
    db[i]=0
    b=b-db
    b[i]=b[i]/u[i,i]
  }
  return=b
}

linearsolves=function(y,X){
  m=as(X,"sparseMatrix")
  m=crossprod(m,m)
  b=crossprod(X,y)
  return=solve(m,b)
}