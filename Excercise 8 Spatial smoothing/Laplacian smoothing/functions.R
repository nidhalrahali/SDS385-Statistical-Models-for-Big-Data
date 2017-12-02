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