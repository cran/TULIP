cv.SeSDA<-function(x, y, nfolds = 5,lambda=NULL, lambda.opt="min", standardize=FALSE, alpha=1, eps=1e-7){
  
  n1<-sum(y==1)
  n2<-sum(y==2)
  
  n<-nrow(x)
  p<-ncol(x)
  x.norm<-matrix(0,n,p)
  transform <- array(list(),p)
  for(i in 1:p){
    obj.norm<-getnorm(x[,i],y)
    x.norm[,i]<-obj.norm$x.norm
    transform[[i]]=obj.norm$transform
  }
  
  
  obj<-cv.dsda(x.norm, y=y, nfolds=nfolds, lambda.opt=lambda.opt, standardize = standardize, lambda = lambda, alpha = alpha, eps = eps)
  
  outlist <- list(transform=transform, objdsda = obj)
  class(outlist) <- c('cv.SeSDA')
  outlist
  
}