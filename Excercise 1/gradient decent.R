grad=function(beta,X,y){
  omega=1/(exp(-X%*%beta)+1)
  return=omega*(1-omega)*(2*y-1)/((2*y-1)*omega+1-y)%*%X
}

gradientdecent=function(X,y,beta0,eps){
  beta=beta0
  for(i in 1:ite){
    beta=beta+eps*grad(beta,X,y)
  }
  return=beta
}