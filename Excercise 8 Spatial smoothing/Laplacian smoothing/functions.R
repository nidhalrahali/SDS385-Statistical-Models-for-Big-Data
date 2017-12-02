library(Matrix)

makeC_sparse = function (dim1, dim2,lambda)  {
  n=dim1*dim2
  D1 = bandSparse(n,  k = c(0, 1), 
                  diagonals = list(rep(-1, n), rep(1, n - 1)))
  D1 = D1[(seq(1, n)%%dim1) != 0, ]
  D2 = bandSparse(n - dim1, m = n, k = c(0, dim1), diagonals = list(rep(-1, n), rep(1, n - 1)))
  D=rBind(D1,D2)
  crossprod(D,D)*lambda+bandSparse(n,k=0,diagonals=list(rep(1,n)))
}

gauss_seidel=function(C,b,ite,targetfunction=FALSE){
  L=tril(C)
  U=C-L
  n=length(b)
  x=solve(L,b,sparse=TRUE)
  tf=rep(0,ite)
  for (i in 1:ite){
    x=solve(L,b-U%*%x,sparse=TRUE)
    if(targetfunction==TRUE)tf[i]=target_function(C,y,x)
  }
  list(x,tf)
}
jacobi=function(C,b,ite,targetfunction=FALSE){
  R=C-band(C,0,0)
  D=diag(C)
  n=length(b)
  bD=b/D
  x=bD
  RD=R/D
  tf=rep(0,ite)
  for (i in 1:ite){
    x=bD-RD%*%x
    if(targetfunction==TRUE)tf[i]=target_function(C,y,x)
  }
  list(x,tf)
}
target_function=function(C,y,x){
  crossprod(x,C%*%x)/2-crossprod(y,x)
}