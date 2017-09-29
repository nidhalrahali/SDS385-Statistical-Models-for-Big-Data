library(Matrix)
linearsolve=function(y,X){
  m=crossprod(X,X)
  b=crossprod(X,y)
  return=solve(m)%*%b
}

linearsolvelu=function(y,X){
  m=crossprod(X,X)
  b=crossprod(X,y)
  lufact=expand(lu(m))
  b=lufact$P%*%b
  b=solve(lufact$L,b)
  b=solve(lufact$U,b)
  return=b
}

linearsolvesparse=function(y,X){
  m=as(X,"sparseMatrix")
  m=crossprod(m,m)
  b=crossprod(X,y)
  return=solve(m,b,sparse=TRUE)
}