SOS<-function(x, y, standardize = FALSE, lambda = NULL, eps = 1e-07){
  obj <- dsda(x, y=y, standardize=standardize, lambda=lambda, alpha=1, eps=eps)
  p=dim(x)[2]
  n=length(y)
  nlambda=length(obj$lambda)
  beta=obj$beta[2:(p+1),]
  lambda=obj$lambda
  
  
  pi1=sum(y==1)/n
  pi2=sum(y==2)/n
  newlambda=lambda*sqrt(pi1*pi2)
  newbeta=beta*sqrt(pi1*pi2)
  
  outlist=list(beta=newbeta, lambda=newlambda)
  class(outlist)=c("SOS")
  outlist
}