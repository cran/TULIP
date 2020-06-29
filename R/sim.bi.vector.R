sim.bi.vector <- function(tesize = 100){
#generate training data
#no kronecker product structure
  p <- 500
  k <- 2
  rho1 <- 0.7
  rho3 <- 0.3

  sigma <- matrix(0,p,p)

  for (i in 1:p){
    for (j in 1:p){
      if (i!=j)
        sigma[i,j] <- rho3
      else
        sigma[i,i] <- 1
    }
  }

  
  nk <- 75
  n <- nk*k

  #define dsigma 
  dsigma <- t(chol(sigma))


  #define B and mean
  beta <- matrix(0,nrow=p,ncol=1)
  beta[1:10,1] <- 0.5
  M <- matrix(0,nrow=2,ncol=p)
  M[2,] <- sigma%*%beta


  y <- rep(0,n)
  for (i in 1:k){
    y[(nk*(i-1)+1):(nk*i)] <- i
  }


  #generate test data

  telabel <- ceiling(runif(tesize)*k)

  vec_x <- matrix(rnorm(p*n),ncol=p)
  x <- vec_x%*%t(dsigma)
  x[y==2,] <- x[y==2,]+M[2,]

  vec_testx <- matrix(rnorm(p*tesize),ncol=p)
  testx <- vec_testx%*%t(dsigma)
  testx[telabel==2,] <- testx[telabel==2,]+M[2,]

  
  outlist <- list(x=x, testx=testx, y=y, testy=telabel)
}
