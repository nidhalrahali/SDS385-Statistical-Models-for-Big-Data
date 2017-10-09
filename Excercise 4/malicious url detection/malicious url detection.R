library(readr)
library(Matrix)
library(microbenchmark)
source('~/GitHub/SDS385-course-work/Excercise 4/malicious url detection/gradient decent functions.R')
source('~/GitHub/SDS385-course-work/Excercise 4/malicious url detection/gradient decent functions cpp.R')
Xtest=readRDS('~/url_svmlight/url_Xtest.rds')
ytest=readRDS('~/url_svmlight/url_ytest.rds')
Xtrain=readRDS('~/url_svmlight/url_Xtrain.rds')
ytrain=readRDS('~/url_svmlight/url_ytrain.rds')

ite = 1
eps = 1.0
lambda = 1.0
beta0 = rep(0,Xtrain@Dim[2])
rn = sample(length(ytrain),ite)

t1=microbenchmark(sgd_adagrad(rn,Xtrain,ytrain,beta0,eps,ite,lambda),times=1L)
result = sgd_adagrad(rn,Xtrain,ytrain,beta0,eps,ite,lambda)
plot(result[[2]],type='l',xlab='',ylab='negative loglikelihood')
ogtest = omega(Xtest,result[[1]])
nllh(ogtest,ytest)/length(ytest)

t2=microbenchmark( sgdC(rn, Xtrain ,ytrain,beta0,eps,ite,lambda),times=1L)
result_cpp = sgdC(rn, Xtrain ,ytrain,beta0,eps,ite,lambda)
plot(result[1:ite],type='l')
beta = result[(ite+1):(ite+Xtrain@Dim[2])]
ogtest = omega(Xtest,beta)
nllh(ogtest,ytest)/length(ytest)
