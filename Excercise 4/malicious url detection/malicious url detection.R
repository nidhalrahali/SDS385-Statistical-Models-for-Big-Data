library(readr)
library(Matrix)
source('~/GitHub/SDS385-course-work/Excercise 4/malicious url detection/gradient decent functions.R')
source('~/GitHub/SDS385-course-work/Excercise 4/malicious url detection/gradient decent functions cpp.R')
Xtest=readRDS('~/url_Xtest.rds')
ytest=readRDS('~/url_ytest.rds')
Xtrain=readRDS('~/url_Xtrain.rds')
ytrain=readRDS('~/url_ytrain.rds')

ite=10
eps=1.0
lambda=1.0
beta0=rep(0,Xtrain@Dim[2])

result=sgd_adagrad(Xtrain,ytrain,beta0,eps,ite,lambda)

ogtest=omega(Xtest,result[[1]])
nllh(ogtest,ytest)/length(ytest)

rn=sample(length(ytrain),ite)
result = sgdC(rn, Xtrain ,ytrain,beta0,eps,ite,lambda)
plot(result[1:ite],type='l')
beta = result[(ite+1):(ite+Xtrain@Dim[2])]
ogtest = omega(Xtest,beta)
nllh(ogtest,ytest)/length(ytest)
