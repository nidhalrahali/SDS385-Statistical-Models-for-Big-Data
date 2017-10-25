library(Matrix)

# solve linear system using matrix inversion
linearsolve=function(y,X){
  m=crossprod(X,X)
  b=crossprod(X,y)
  return=solve(m)%*%b
}

# solves linear system by first using LU decomposition, and then solving two linear systems with
# optimized method for triangular matrix
linearsolvelu=function(y,X){
  m=crossprod(X,X)
  b=crossprod(X,y)
  lufact=expand(lu(m))
  b=lufact$P%*%b
  b=solve(lufact$L,b)
  b=solve(lufact$U,b)
  return=b
}

# solve linear system by optimized method for sparse matrix
linearsolvesparse=function(y,X){
  m=as(X,"sparseMatrix")
  m=crossprod(m,m)
  b=crossprod(X,y)
  return=solve(m,b,sparse=TRUE)
}