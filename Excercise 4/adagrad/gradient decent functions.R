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

sgd=function(X,y,beta0,eps,ite){
  beta=beta0
  betahistory=matrix(nrow=length(beta0),ncol=ite)
  for(i in 1:ite){
    r = sample(length(y),1)
    og = omega(X[r,],beta)
    beta = beta-eps*grad(og,y[r],X[r,])
    betahistory[,i] = beta
  }
  betahistory
}

sgd_adagrad=function(X,y,beta0,eps,ite){
  beta=beta0
  H=rep(0,length(beta))
  betahistory=matrix(nrow=length(beta0),ncol=ite)
  for(i in 1:ite){
    r = sample(length(y),1)
    og = omega(X[r,],beta)
    g=grad(og,y[r],X[r,])
    H=H+g^2
    beta = beta-eps*g/sqrt(H)
    betahistory[,i] = beta
  }
  betahistory
}
