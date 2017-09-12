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

newtondirection=function(omega,y,X){
  w=omega*(1-omega)*(2*y-1)*((1-2*omega)*((2*y-1)*omega+1-y)-(2*y-1)*omega*(1-omega))/((2*y-1)*omega+1-y)^2
  XW=X
  for(i in 1:nrow(X)){
    XW[i,]=XW[i,]*w[i]
  }
  g=grad(omega,y,X)
  return=-solve(crossprod(X,XW))%*%g
}

newtonmethod=function(X,y,beta0,ite){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    omega=as.vector(1/(exp(-X%*%beta)+1))
    nllh[i]=nllh(omega,y)
    beta=beta+newtondirection(omega,y,X)
  }
  return=list(beta=beta,negloglikelihood=nllh)
}