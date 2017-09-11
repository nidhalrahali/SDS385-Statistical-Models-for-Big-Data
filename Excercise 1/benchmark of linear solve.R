library(Matrix)

source('~/GitHub/SDS385-course-work/Excercise 1/linear solver.R')

n=5000
p=2000

X=matrix(nrow=n,ncol=p)
beta=matrix(1,nrow=p,ncol=1)
for(i in 1:n){
  for(j in 1:p)X[i,j]=runif(1,min=0,max=1)
}
y=X%*%beta+rnorm(1,sd=10)
W=matrix(1,nrow=n)

betas=linearsolve1(y,X,W)

m=crossprod(X,X)
lufact=expand(lu(m))
l=lufact$L

system.time(linearsolve1(y,X,W))
betas=linearsolve2(y,X,W)
