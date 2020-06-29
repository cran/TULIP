sim.tensor.cov <- function(tesize = 100){
#generate training data
#no kronecker product structure
  dimen <- c(10,10,10)
  nvars <- prod(dimen)
  k <- 2
  rho1 <- 0.7
  rho3 <- 0.3

  sigma <- array(list(),length(dimen))
  dsigma <- array(list(),length(dimen))

  sigma[[1]] <- diag(dimen[1])
  sigma[[2]] <- diag(dimen[2])
  sigma[[3]] <- diag(dimen[3])

  nk <- 75
  n <- nk*k

  #define dsigma 
  for (i in 1:length(dimen)){
    dsigma[[i]] <- t(chol(sigma[[i]]))
  }


  #define B and mean
  B2 <- array(0,dim=dimen)
  ini <- c(1,2)
  inj <- c(1,2)
  ins <- c(1,2)
  for (i in ini){
    for (j in inj){
      for (s in ins){
        B2[i,j,s] <- 0.8
      }
    }
  }
  M <- array(list(),k)
  M[[1]] <- array(0,dim=dimen)
  M[[2]] <- atrans(B2,sigma)


  y <- rep(0,n)
  for (i in 1:k){
    y[(nk*(i-1)+1):(nk*i)] <- i
  }


  coef <- array(0,dim=c(dimen,2))
  for(i in 1:5){
	  for (j in 1:5){
		  for (s in 1:5){
			  coef[i,j,s,1] <- 1
		  }
	  }
  }

  #generate test data

  telabel <- ceiling(runif(tesize)*k)

  z <- matrix(rnorm(2*n),nrow=n,ncol=2)
  z[y==1,] <- z[y==1,]
  z[y==2,] <- z[y==2,]+0.3


  testz <- matrix(rnorm(2*tesize),nrow=tesize,ncol=2)
  testz[telabel==1,] <- testz[telabel==1,]
  testz[telabel==2,] <- testz[telabel==2,]+0.3

  vec_x <- matrix(rnorm(nvars*n),ncol=n)
  x <- array(list(),n)
  for (i in 1:n){
      x[[i]] <- array(vec_x[,i],dimen)
      x[[i]] <- x[[i]]+amprod(coef,t(z[i,]),4)[,,,1]
      x[[i]] <- M[[y[i]]]+atrans(x[[i]],dsigma)
  }


  vec_testx <- matrix(rnorm(nvars*tesize),ncol=tesize)
  testx <- array(list(),tesize)
  for (i in 1:tesize){
    testx[[i]] <- array(vec_testx[,i],dimen)
    testx[[i]] <- testx[[i]]+amprod(coef,t(testz[i,]),4)[,,,1]
    testx[[i]] <- M[[telabel[i]]]+atrans(testx[[i]],dsigma)
  }
  
  vec_x <- t(vec_x)
  vec_testx <- t(vec_testx)
  
  
  outlist <- list(x=x, z=z, testx=testx, testz=testz, vec_x=vec_x, vec_testx=vec_testx, y=y, testy=telabel)
}
