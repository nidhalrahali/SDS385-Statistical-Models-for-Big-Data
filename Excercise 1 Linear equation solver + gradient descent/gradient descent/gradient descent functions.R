# given observation X and beta, compute logistic function omega

omega=function(X,beta){
  return=as.vector(1/(exp(-X%*%beta)+1))
}

# given omega, observation y and X,

grad=function(omega,y,X){
  return=as.vector(crossprod(omega-y,X))
}

# compute negative loglikelihood given omega and y
nllh=function(omega,y){
  r=0
  for(i in 1:length(y)){
    if(y[i]==1)r=r-log(omega[i])
    else r=r-log(1-omega[i])
  }
  return=r
}
#run gradient decent with training data X y and test data testX testy,step size eps and iteration number ite
gradientdecent=function(X,y,testX,testy,beta0,eps,ite){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  tnllh=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    og=omega(X,beta)
    testog=omega(testX,beta)
    nllh[i]=nllh(og,y)/length(y)
    tnllh[i]=nllh(testog,testy)/length(testy)
    beta=beta-eps*grad(og,y,X)
  }
  return=list(beta=beta,negloglikelihood=nllh,testnegloglikelihood=tnllh)
}

# compute newton direction in the space of beta
newtondirection=function(omega,y,X){
  w=omega*(1-omega)
  g=grad(omega,y,X)
  return=as.vector(-solve(crossprod(X,X*w),g))
}

# run newton direction
newtonmethod=function(X,y,testX,testy,beta0,ite){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  tnllh=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    og=omega(X,beta)
    testog=omega(testX,beta)
    nllh[i]=nllh(og,y)/length(y)
    tnllh[i]=nllh(testog,testy)/length(testy)
    beta=beta+newtondirection(og,y,X)
  }
  return=list(beta=beta,negloglikelihood=nllh,testnegloglikelihood=tnllh)
}