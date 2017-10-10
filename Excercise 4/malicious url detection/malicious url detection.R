library(readr)
library(Matrix)
library(microbenchmark)
library(Rcpp)
library(RcppEigen)
Rcpp::sourceCpp('gradient_decent_functions.cpp')
source('~/GitHub/SDS385-course-work/Excercise 4/malicious url detection/gradient decent functions.R')
Xtest=readRDS('~/url_svmlight/url_Xtest.rds')
ytest=readRDS('~/url_svmlight/url_ytest.rds')
Xtrain=readRDS('~/url_svmlight/url_Xtrain.rds')
Xtrain=t(Xtrain)
ytrain=readRDS('~/url_svmlight/url_ytrain.rds')

eps = 1
lambda = 1
rn = sample(length(ytrain))
epoch=1

t1=microbenchmark( sgdC(rn, Xtrain ,ytrain,eps,ite,lambda),times=1L)
result = sgdC(rn, Xtrain ,ytrain,eps,ite,lambda)
plot(result[1:ite],type='l')
beta = result[(ite+1):(ite+Xtrain@Dim[2])]
ogtest = omega(Xtest,beta)
nllh(ogtest,ytest)/length(ytest)

t2=microbenchmark(sgdC_sparse(rn,Xtrain,ytrain,eps,ite,lambda),times=1L)
result_sparse=sgdC_sparse(rn,Xtrain,ytrain,eps,epoch,lambda)
plot(result_sparse[1:ncol(Xtrain)*epoch],type='l')
