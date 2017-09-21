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
beta0=as.vector(matrix(0,nrow=11))

trainX=X[1:250,]
trainy=y[1:250]
testX=X[251:569,]
testy=y[251:569]

ite=1000
eps=0.02
alpha=0.05
result1=stochasticgradientdecent(trainX,trainy,testX,testy,beta0,eps,ite,alpha)
plot(result1$averagenegloglikelihood,type='l')

eps=0.05
result2=stochasticgradientdecent(trainX,trainy,testX,testy,beta0,eps,ite,alpha)
plot(result2$averagenegloglikelihood,type='l')

eps=0.1
result3=stochasticgradientdecent(trainX,trainy,testX,testy,beta0,eps,ite,alpha)
plot(result3$averagenegloglikelihood,type='l')
par(mfrow=c(2,2))
plot(result1$negloglikelihood,type='l',ylim=c(0.1,0.7))
lines(result1$testnegloglikelihoood,col='red')

plot(result2$negloglikelihood,type='l',ylim=c(0.1,0.7))
lines(result2$testnegloglikelihoood,col='red')

plot(result3$negloglikelihood,type='l',ylim=c(0.1,0.7))
lines(result3$testnegloglikelihoood,col='red')