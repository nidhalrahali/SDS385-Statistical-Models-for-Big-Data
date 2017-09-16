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
beta0=as.vector(matrix(1,nrow=11))

trainX=X[1:250,]
trainy=y[1:250]

ite=1000
eps=0.01
result1=gradientdecent(trainX,trainy,beta0,eps,ite)
plot(result1$negloglikelihood)
result1$negloglikelihood
result1$beta

testX=X[251:569,]
testy=y[251:569]
likelihood=-nllh(1/(exp(-as.vector(testX%*%result1$beta))+1),testy)


ite=20
beta0=as.vector(matrix(0.7,nrow=11))
result2=newtonmethod(trainX,trainy,beta0,ite)
result2$beta
plot(result2$negloglikelihood)
result2$negloglikelihood
likelihood=-nllh(1/(exp(-as.vector(testX%*%result2$beta))+1),testy)

result3=gradientdecent_varyingstep(X,y,beta0,ite,0.7,0.9)
plot(result3$negloglikelihood)
result3$negloglikelihood
result3$beta
result3$step