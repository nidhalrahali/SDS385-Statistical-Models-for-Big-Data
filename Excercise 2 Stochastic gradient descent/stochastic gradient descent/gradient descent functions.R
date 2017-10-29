omega=function(X,beta){
  return=as.vector(1/(exp(-X%*%beta)+1))
}

grad=function(omega,y,X){
  return=as.vector(crossprod(omega-y,X))
}

nllh=function(omega,y){
  r=0
  for(i in 1:length(y)){
    if(y[i]==1)r=r-log(omega[i])
    else r=r-log(1-omega[i])
  }
  return=r
}

#run stochastic gradient decent on training data X y and test data testX,testY using with
# iterations number=ite and step size=eps
# evaluate the moving average of negative loglikelihood with exponential decay rate alpha
stochasticgradientdescent=function(X,y,testX,testy,beta0,eps,ite,alpha){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  tnllh=as.vector(matrix(nrow=ite))
  nllhaverage=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    # draw the index of direction
    r=sample(length(y),1)
    og=omega(X,beta)
    testog=omega(testX,beta)
    # we still compute total negative loglikelihood at each step for reference, this is not
    # necessary in practice
    nllh[i]=nllh(og,y)/length(y)
    tnllh[i]=nllh(testog,testy)/length(testy)
    beta=beta-eps*grad(og[r],y[r],X[r,])
    if(i==1)nllhaverage[i]=nllh(og[r],y[r])
    else nllhaverage[i]=nllh(og[r],y[r])*alpha+nllhaverage[i-1]*(1-alpha)
  }
  return=list(beta=beta,negloglikelihood=nllh,testnegloglikelihoood=tnllh,averagenegloglikelihood=nllhaverage)
}
#run stochastic gradient decent with vary steps, the steps are computed using Robbins Monro rule
# step size at t=C/(t0+t)^decay
varyingstepsgradientdescent=function(X,y,testX,testy,beta0,ite,alpha,decay,t0,C){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  tnllh=as.vector(matrix(nrow=ite))
  nllhaverage=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    r=sample(length(y),1)
    og=omega(X,beta)
    testog=omega(testX,beta)
    nllh[i]=nllh(og,y)/length(y)
    tnllh[i]=nllh(testog,testy)/length(testy)
    eps=C/((t0+i)^decay)
    beta=beta-eps*grad(og[r],y[r],X[r,])
    if(i==1)nllhaverage[i]=nllh(og[r],y[r])
    else nllhaverage[i]=nllh(og[r],y[r])*alpha+nllhaverage[i-1]*(1-alpha)
  }
  return=list(beta=beta,negloglikelihood=nllh,testnegloglikelihoood=tnllh,averagenegloglikelihood=nllhaverage)
}