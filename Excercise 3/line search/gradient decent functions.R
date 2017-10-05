# given observation X and beta, compute logistic function omega

backtrack=function(X,y,beta,og,direction,eps0,alpha,c){
  current_nllh=nllh(og,y)
  sufficient=FALSE
  eps=eps0
  while(!sufficient){
    eps=eps*alpha
    propose_beta=beta+eps*direction
    propose_og=omega(X,propose_beta)
    if(nllh(propose_og,y)<current_nllh-c*eps*crossprod(direction,direction))sufficient=TRUE
  }
  eps
}

omega=function(X,beta){
  as.vector(1/(exp(-X%*%beta)+1))
}

# given omega, observation y and X,

grad=function(omega,y,X){
  as.vector(crossprod(omega-y,X))
}

# compute negative loglikelihood given omega and y
nllh=function(omega,y){
  r=0
  for(i in 1:length(y)){
    if(y[i]==1)r=r-log(omega[i])
    else r=r-log(1-omega[i])
  }
  r
}
#run gradient decent with training data X y and test data testX testy,step size eps and iteration number ite
gradientdecent=function(X,y,beta0,eps,ite){
  betahistory=matrix(nrow=length(beta0),ncol=ite)
  beta=beta0
  for(i in 1:ite){
    og = omega(X,beta)
    beta = beta-eps*grad(og,y,X)
    betahistory[,i] = beta
  }
  betahistory
}

gradientdecent_linesearch=function(X,y,beta0,eps0,ite,alpha,c){
  betahistory = matrix(nrow=length(beta0),ncol=ite)
  beta = beta0
  epshistory = rep(0,ite)
  for(i in 1:ite){
    og = omega(X,beta)
    direction=-grad(og,y,X)
    eps=backtrack(X,y,beta,og,direction,eps0,alpha,c)
    epshistory[i]=eps
    beta = beta+eps*direction
    betahistory[,i] = beta
  }
  list(betahistory=betahistory,epshistory=epshistory)
}
# compute newton direction in the space of beta
newtondirection=function(omega,y,X){
  w=omega*(1-omega)
  g=grad(omega,y,X)
  as.vector(-solve(crossprod(X,X*w),g))
}

# run newton direction
newtonmethod=function(X,y,beta0,ite){
  beta=beta0
  betahistory = matrix(nrow = length(beta0),ncol = ite)
  for(i in 1:ite){
    og=omega(X,beta)
    beta=beta+newtondirection(og,y,X)
    betahistory[,i] = beta
  }
  betahistory=betahistory
}

quasi_newtonmethod=function(X,y,beta0,eps0,ite,alpha,c){
  beta = beta0
  H = diag(1,nrow = length(beta0),ncol = length(beta0))
  betahistory = matrix(nrow = length(beta0),ncol = ite)
  for(i in 1:ite){
    if(i==1){
      og = omega(X,beta)
      g = grad(og,y,X)
    }
    else{
      og=ognew
      g=gnew
    }
    direction=as.vector(-H%*%g)
    eps=backtrack(X,y,beta,og,direction,eps0,alpha,c)
    e = eps*direction
    beta = beta+e
    betahistory[,i] = beta
    ognew = omega(X,beta)
    gnew = grad(ognew,y,X)
    if(crossprod(gnew,gnew)<1)break
    t = gnew-g
    Ht=H%*%t
    H = H- tcrossprod(Ht,Ht)/as.numeric(crossprod(Ht,t)) + tcrossprod(e,e)/as.numeric(crossprod(t,e))
  }
  list(betahistory=betahistory)
}