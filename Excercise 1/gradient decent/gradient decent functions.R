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

newtondirection=function(omega,y,X){
  w=omega*(1-omega)
  g=grad(omega,y,X)
  return=as.vector(-solve(crossprod(X,X*w),g))
}

newtonmethod=function(X,y,beta0,ite){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    og=omega(X,beta)
    nllh[i]=nllh(og,y)
    beta=beta+newtondirection(og,y,X)
  }
  return=list(beta=beta,negloglikelihood=nllh)
}