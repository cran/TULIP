SeSDA<-function(x, y, standardize = FALSE, lambda = NULL, alpha = 1, eps = 1e-07){
  
  n<-nrow(x)
  p<-ncol(x)
  x.norm<-matrix(0,n,p)
  transform <- array(list(),p)
  for(i in 1:p){
    obj.norm<-getnorm(x[,i],y)
    x.norm[,i]<-obj.norm$x.norm
    transform[[i]]=obj.norm$transform
  }
  
  obj<-dsda(x.norm, y=y, standardize = standardize, lambda = lambda, alpha = alpha, eps = eps)
  
  outlist <- list(transform=transform, objdsda = obj)
  class(outlist) <- c('SeSDA')
  outlist
  
}