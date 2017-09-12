grad=function(omega,y,X){
  return=as.vector(crossprod(omega*(1-omega)*(2*y-1)/((2*y-1)*omega+1-y),X))
}

nllh=function(omega,y){
  return=-sum(log((2*y-1)*omega+1-y))
}

gradientdecent=function(X,y,beta0,eps,ite){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    omega=as.vector(1/(exp(-X%*%beta)+1))
    nllh[i]=nllh(omega,y)
    beta=beta+eps*grad(omega,y,X)
  }
  return=list(beta=beta,negloglikelihood=nllh)
}