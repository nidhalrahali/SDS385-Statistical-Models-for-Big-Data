target_function=function(X,Y,beta,lambda){
  mean((Y-X%*%beta)^2)+lambda*sum(abs(beta[2:length(beta)]))
}

soft_threshold=function(y,lambda){
  ret=as.vector(y)
  for(i in 2:length(y)){
    ret[i] = sign(y[i])*max(c(0,abs(y[i])-lambda))
  }
  ret
}

gradient=function(X,Y,beta){
  crossprod(X,X%*%beta-Y)/length(Y)
}


acc_prox_gd=function(X,Y,lambda,gamma,ite){
  tf_history=rep(0,ite)
  p=ncol(X)
  s1=10
  z=rep(0,p)
  beta1=rep(0,p)
  for(it in 1:ite){
    grad=gradient(X,Y,z)
    u=z-gamma*grad
    beta2=soft_threshold(u,lambda*gamma)
    tf_history[it]=target_function(X,Y,beta2,lambda)
    s2=(1+sqrt(1+4*s1^2))/2
    z=beta2+(s1-2)/s2*(beta2-beta1)
    beta1=beta2
    s1=s2
  }
  list(beta=beta2,targetfunction=tf_history)
}

admm_lasso=function(X,Y,lambda,rho,ite){
  tf_history=rep(0,ite)
  n=length(Y)
  p=ncol(X)
  u=rep(0,p)
  beta=rep(0,p)
  alpha=rep(0,p)
  XX=solve(crossprod(X,X)/n+rho*diag(nrow=p,ncol=p))
  XY=as.vector(crossprod(X,Y)/n)
  lambda_over_rho=lambda/rho
  for(it in 1:ite){
    tf_history[it]=target_function(X,Y,beta,lambda)
    v=alpha-u
    beta=as.vector(XX%*%(XY+rho*v))
    w=beta+u
    alpha=soft_threshold(w,lambda_over_rho)
    u=u+beta-alpha
  }
  list(beta=beta,targetfunction=tf_history)
}