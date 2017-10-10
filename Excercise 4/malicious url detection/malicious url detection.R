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

t=microbenchmark(sgdC_sparse(rn,Xtrain,ytrain,eps,epoch,lambda),times=1L)
result=sgdC_sparse(rn,Xtrain,ytrain,eps,epoch,lambda)
plot(result[1:ncol(Xtrain)*epoch],type='l')
beta=result[(ncol(Xtrain)*epoch+1):length(result)]

#prediction
og=omega(Xtest,beta)
ypredict={og>0.5}
dif=ypredict-ytest
#true positive
tp=sum({dif==0&&ypredit==1})/length(ytest)
#true negative
tp=sum({dif==0&&ypredit==0})/length(ytest)
#false positive
fp=sum({dif==1})/length(ytest)
#false negative
fn=sum({dif==-1})/length(ytest)