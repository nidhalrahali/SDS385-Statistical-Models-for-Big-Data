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

stochasticgradientdecent=function(X,y,testX,testy,beta0,eps,ite,alpha){
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
    beta=beta-eps*grad(og[r],y[r],X[r,])
    if(i==1)nllhaverage[i]=nllh(og[r],y[r])
    else nllhaverage[i]=nllh(og[r],y[r])*alpha+nllhaverage[i-1]*(1-alpha)
  }
  return=list(beta=beta,negloglikelihood=nllh,testnegloglikelihoood=tnllh,averagenegloglikelihood=nllhaverage)
}

gradientdecent=function(X,y,beta0,eps,ite){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    og=omega(X,beta)
    nllh[i]=nllh(og,y)
    beta=beta-eps*grad(og,y,X)
  }
  return=list(beta=beta,negloglikelihood=nllh)
}