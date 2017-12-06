library(Matrix)

makeD_sparse = function (dim1, dim2)  {
  n=dim1*dim2
  D1 = bandSparse(n,  k = c(0, 1), 
                  diagonals = list(rep(-1, n), rep(1, n - 1)))
  D1 = D1[(seq(1, n)%%dim1) != 0, ]
  D2 = bandSparse(n - dim1, m = n, k = c(0, dim1), diagonals = list(rep(-1, n), rep(1, n - 1)))
  rBind(D1,D2)
}

soft_threshold=function(y,lambda){
  sign(y)*max(c(0,abs(y)-lambda))
}

admm_lasso=function(y,D,lambda,rho,ite){
  tf_history=rep(0,ite)
  n=ncol(D)
  e=nrow(D)
  u=rep(0,e)
  x=rep(0,n)
  z=rep(0,e)
  C=crossprod(D,D)*lambda+bandSparse(n,k=0,diagonals=list(rep(1,n)))
  lufact=expand(lu(C))
  lambda_over_rho=lambda/(2*rho)
  Dx=rep(0,e)
  for(it in 1:ite){
    tf_history[it]=sum((y-x)^2)/2+lambda*sum(abs(Dx))/2
    v=z-u
    x=solve(lufact$L,lufact$P%*%(y+rho*crossprod(D,z)),sparse=TRUE)
    x=solve(lufact$Q)%*%solve(lufact$U,x,sparse=TRUE)
    Dx=D%*%x
    w=Dx+u
    z=apply(w,1,soft_threshold,lambda=lambda_over_rho)
    u=u+Dx-z
  }
  list(x=x,targetfunction=tf_history)
}