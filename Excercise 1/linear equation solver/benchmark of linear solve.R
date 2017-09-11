library(Matrix)
library(microbenchmark)
source('~/GitHub/SDS385-course-work/Excercise 1/linear solver.R')

n=5000
p=2000

X=matrix(runif(n*p,0,1),nrow=n,ncol=p)
beta=matrix(1,nrow=p,ncol=1)
y=X%*%beta+rnorm(1,sd=10)


l1=microbenchmark(linearsolve1(y,X),times=1L)
l2=microbenchmark(linearsolve2(y,X),times=1L)

mask=matrix(rbinom(n*p,1,0.05),nrow=n,ncol=p)
Xsparse=X*mask
y=Xsparse%*%beta+rnorm(1,sd=10)
ls=microbenchmark(linearsolves(y,Xsparse),times=1L)
l1=microbenchmark(linearsolve1(y,Xsparse),times=1L)
l2=microbenchmark(linearsolve2(y,Xsparse),times=1L)


