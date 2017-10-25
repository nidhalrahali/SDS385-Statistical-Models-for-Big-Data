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


prox_gd=function(X,Y,beta0,lambda,gamma,ite){
  beta=beta0
  tf_history=rep(0,ite)
  for(it in 1:ite){
    grad=gradient(X,Y,beta)
    u=beta-gamma*grad
    beta=soft_threshold(u,lambda*gamma)
    tf_history[it]=target_function(X,Y,beta,lambda)
  }
  list(beta=beta,targetfunction=tf_history)
}