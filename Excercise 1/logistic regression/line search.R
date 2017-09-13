grad=function(omega,y,X){
  return=as.vector(crossprod(omega*(1-omega)*(2*y-1)/((2*y-1)*omega+1-y),X))
}

nllh=function(omega,y){
  return=-sum(log((2*y-1)*omega+1-y))
}

wolfcondition=function(eps, beta, g,y,X,c1,c2){
  omega=as.vector(1/(exp(-X%*%beta)+1))
  nllhoriginal=nllh(omega,y)
  betatrial=beta+eps*g
  omega=as.vector(1/(exp(-X%*%betatrial)+1))
  nllhtrial=nllh(omega,y)
  gtrial=grad(omega,y,X)
  if(nllhtrial<nllhoriginal+c1*eps*crossprod(g,beta) & abs(crossprod(gtrial,beta))<c2*abs(crossprod(g,beta)))return=TRUE
  else return=FALSE
}

gradientdecent_varyingstep=function(X,y,beta0,ite,c1,c2){
  beta=beta0
  nllh=as.vector(matrix(nrow=ite))
  accepteps=as.vector(matrix(nrow=ite))
  for(i in 1:ite){
    omega=as.vector(1/(exp(-X%*%beta)+1))
    nllh[i]=nllh(omega,y)
    g=grad(omega,y,X)
    for(j in 1:1000){
      eps=j/5000
      if(wolfcondition(eps,beta,g,y,X,c1,c2))break
    }
    accepteps[i]=eps
    beta=beta+eps*g
  }
  return=list(beta=beta,negloglikelihood=nllh,step=accepteps)
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
  return=as.vector(-solve(crossprod(X,XW))%*%g)
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