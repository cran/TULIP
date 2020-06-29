adjvec <- function(x, z, y, testx=NULL, testz=NULL, is.centered=FALSE){
  
  
  n <- length(y)
  z <- as.matrix(z)
  p <- dim(x)[2]
  q <- dim(z)[2]
  
  cz <- z
  cx <- x
  nclass <- as.integer(length(unique(y)))
  if (is.centered==FALSE){
    for (i in 1:nclass){
      if (q>1) {cz[y == i,] <- sweep(z[y == i,],2,colMeans(z[y == i,]))}
      else {cz[y==i] <- z[y==i]-mean(z[y==i])}
      cx[y == i,] <- sweep(x[y == i,],2,colMeans(x[y == i,]))
    }
  }
  c <- solve(t(cz) %*% cz) %*% t(cz)
  c <- c %*% cx
  
  xres <- x-z%*%c
  
  if (!is.null(testx)){
    testxres <- testx-testz %*% c
  }else{
    testxres <- NULL
  }
  
  muz <- matrix(0,nrow=q,ncol=(nclass-1))
  for (i in 2:nclass){
    if (q>1){muz[,(i-1)] <- apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)}
    else {muz[i-1] <- mean(z[y==i])-mean(z[y==1])}
  }
  gamma <- solve(cov(z)) %*% muz
  
  outlist = list(gamma=as.matrix(gamma),xres=xres,testxres=testxres)
  outlist
}