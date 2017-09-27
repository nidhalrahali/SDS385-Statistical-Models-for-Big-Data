library(readr)
data <- read_csv("~/GitHub/SDS385-course-work/Excercise 1/logistic regression/wdbc.csv",col_names = FALSE)
source('~/GitHub/SDS385-course-work/Excercise 2/stochastic logistic regression/gradient decent functions.R')
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

# experiments with different steps, the average negative likelihood is
# exponentially decay moving average
ite=1000
alpha=0.05
eps=0.02
result1=stochasticgradientdecent(trainX,trainy,testX,testy,beta0,eps,ite,alpha)
plot(result1$averagenegloglikelihood,type='l')

eps=0.05
result2=stochasticgradientdecent(trainX,trainy,testX,testy,beta0,eps,ite,alpha)
plot(result2$averagenegloglikelihood,type='l')

eps=0.1
result3=stochasticgradientdecent(trainX,trainy,testX,testy,beta0,eps,ite,alpha)
plot(result3$averagenegloglikelihood,type='l')

#compare the negative likelihood of test data and training data during the iteration
par(mfrow=c(2,2))
plot(result1$negloglikelihood,type='l',ylim=c(0.1,0.7),ylab='negative loglikelihood',xlab='',sub='step=0.02')
lines(result1$testnegloglikelihoood,col='red')

plot(result2$negloglikelihood,type='l',ylim=c(0.1,0.7),ylab='negative loglikelihood',xlab='',sub='step=0.05')
lines(result2$testnegloglikelihoood,col='red')

plot(result3$negloglikelihood,type='l',ylim=c(0.1,0.7),ylab='negative loglikelihood',xlab='',sub='step=0.1')
lines(result3$testnegloglikelihoood,col='red')


#decaying steps
par(mfrow=c(2,2))
t0=1
C=0.5
decay=0.95
result4=varyingstepsgradientdecent(trainX,trainy,testX,testy,beta0,ite,alpha,decay,t0,C)
plot(result4$averagenegloglikelihood,type='l',ylab='average negative loglikelihood',xlab='',sub='alpha=0.95')

decay=0.75
result5=varyingstepsgradientdecent(trainX,trainy,testX,testy,beta0,ite,alpha,decay,t0,C)
plot(result5$averagenegloglikelihood,type='l',ylab='average negative loglikelihood',xlab='',sub='alpha=0.75')

decay=0.6
result6=varyingstepsgradientdecent(trainX,trainy,testX,testy,beta0,ite,alpha,decay,t0,C)
plot(result6$averagenegloglikelihood,type='l',ylab='average negative loglikelihood',xlab='',sub='alpha=0.6')
