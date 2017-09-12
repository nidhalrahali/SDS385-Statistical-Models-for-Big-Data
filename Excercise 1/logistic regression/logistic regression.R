library(readr)
data <- read_csv("~/GitHub/SDS385-course-work/Excercise 1/logistic regression/wdbc.csv",col_names = FALSE)
source('~/GitHub/SDS385-course-work/Excercise 1/logistic regression/line search.R')
X=as.matrix(data[3:12])
for(i in 1:ncol(X)){
  X[,i]=(X[,i]-mean(X[,i]))/sqrt(var(X[,i]))
}
X=cbind(X,1)
y=as.vector(matrix(nrow=nrow(data),ncol=1))
for(i in 1:nrow(data)){
  if(data[i,2]=="M")y[i]=1
  else y[i]=0
}
beta0=as.vector(matrix(0.1,nrow=11))
ite=10000
eps=0.02
resultgd=gradientdecent(X,y,beta0,eps,ite)
plot(resultgd$negloglikelihood)
resultgd$beta

ite=20
resultnt=newtonmethod(X,y,beta0,ite)
resultnt$beta
plot(resultnt$negloglikelihood)
