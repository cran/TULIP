ROAD <- function(x, y, standardize = FALSE, lambda = NULL, eps = 1e-07){
  obj <- dsda(x, y=y, standardize=standardize, lambda=lambda, alpha=1, eps=eps)
  p <- dim(x)[2]
  n <- length(y)
  nlambda <- length(obj$lambda)
  beta <- obj$beta[2:(p+1),]
  lambda <- obj$lambda
  
  newbeta <- beta
  newlambda <- lambda
  
  mu1 <- apply(x[y==1,],2,mean)
  mu2 <- apply(x[y==2,],2,mean)
  w <- mu2-mu1
  for (i in 1:nlambda){
    newlambda[i] <- lambda[i]*2/sum(beta[,i]*w)/n
  }
  beta <- as.matrix(beta)
  newbeta <- as.matrix(newbeta)
  newbeta <- sweep(newbeta, 2, t(beta)%*%w, FUN="/")
  
  outlist <- list(beta=newbeta, lambda=newlambda)
  class(outlist) <- c("ROAD")
  outlist
}