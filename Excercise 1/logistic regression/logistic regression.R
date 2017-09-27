library(readr)
data <- read_csv("~/GitHub/SDS385-course-work/Excercise 1/logistic regression/wdbc.csv",col_names = FALSE)
source('~/GitHub/SDS385-course-work/Excercise 1/logistic regression/gradient decent functions.R')
X=as.matrix(data[3:12])
X=scale(X)
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

ite=200
eps=0.02
result1=gradientdecent(trainX,trainy,beta0,eps,ite)
plot(result1$negloglikelihood,type='l')

#percentage change of likelihood as a measure of convergence
convergence1=abs(result1$negloglikelihood[2:ite]-result1$negloglikelihood[1:(ite-1)])/(result1$negloglikelihood[1:(ite-1)]+0.0001)
plot(convergence1,type='l')

#likelihood for test data
-nllh(1/(exp(-as.vector(testX%*%result1$beta))+1),testy)

#newton method
ite=20
result2=newtonmethod(trainX,trainy,beta0,ite)
plot(result2$negloglikelihood,type='l')
convergence2=abs(result2$negloglikelihood[2:ite]-result2$negloglikelihood[1:(ite-1)])/(result2$negloglikelihood[1:(ite-1)]+0.0001)
plot(convergence2,convergence1)
plot(convergence2[4:ite],convergence1[4:ite],xlab ="newton",ylab="gradient decent" )
#likelihood for training data
-nllh(1/(exp(-as.vector(testX%*%result2$beta))+1),testy)

# a comparison of gradient decent with different step size
ite=200
eps=0.05
result3=gradientdecent(trainX,trainy,beta0,eps,ite)
plot(result3$negloglikelihood,type='l')
convergence3=abs(result3$negloglikelihood[2:ite]-result3$negloglikelihood[1:(ite-1)])/(result3$negloglikelihood[1:(ite-1)]+0.0001)
plot(convergence1[10:ite],convergence3[10:ite],xlab ="eps=0.02",ylab="eps=0.05",xlim=c(0,0.05),ylim=c(0,0.05))
plot(convergence1[50:ite],convergence3[50:ite],xlab ="eps=0.02",ylab="eps=0.05",xlim=c(0,0.0004),ylim=c(0,0.0004))
plot(convergence1[10:ite],type='l')
lines(convergence2[10:ite],col='red')
lines(convergence3[10:ite],col='blue')
