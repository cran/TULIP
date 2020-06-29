predict.msda.original <- function(object, newx, z=NULL, ztest=NULL, gamma=NULL, ...) {

    #newx is the adjusted x if z exists
    beta <- object$beta
    mu <- object$mu
    prior <- object$prior
    n <- nrow(newx)
    p <- ncol(newx)
    x.train <- object$x #object$x is the adjusted x if z exists
    y.train <- object$y
    nclass <- length(prior)
    nlambda <- length(beta)
    pred <- matrix(0, n, nlambda)
    if (!is.null(gamma)){
      betanew <- array(list(),length(beta))
      q <- length(gamma)
      for (i in 1:nlambda){
        betanew[[i]] <- rbind(matrix(beta[[i]],ncol=1),matrix(gamma,ncol=1))
      }
      newxtrain <- cbind(x.train, z)
      newxtest <- cbind(newx, ztest)

      for (i in 1:nlambda){
        xfit <- newxtrain%*%betanew[[i]][,1:(nclass-1),drop=FALSE]
        xfit.sd <- matrix(0, nclass, ncol(xfit))
        for(j in 1:nclass){
          xfit.sd[j,] <- apply(xfit[y.train==j,,drop=FALSE],2,sd)
        }
        xfit.sd <- apply(xfit.sd,2,min)
        if(min(xfit.sd) < 1e-4){
          pred[,i] <- which.max(prior)
        }else{
          l <- lda(xfit, y.train)
          pred[, i] <- predict(l, newxtest %*% betanew[[i]][, 1:(nclass - 1)])$class
        }
      }
    }else{
      for (i in 1:nlambda) {
          nz <- sum(beta[[i]][, 1] != 0)
          if (nz == 0) {
              pred[,i] <- which.max(prior)
          } else {
              xfit <- x.train %*% beta[[i]][, 1:(min(nclass - 1, nz)),drop=FALSE]
              xfit.sd <- matrix(0,nclass,ncol(xfit))
              for(j in 1:nclass){
                xfit.sd[j,] <- apply(xfit[y.train==j,,drop=FALSE],2,sd)
              }
              xfit.sd <- apply(xfit.sd,2,min)
              if(min(xfit.sd) < 1e-4){
                  pred[,i] <- which.max(prior)
                }else{
                l <- lda(xfit, y.train)
                pred[, i] <- predict(l, newx %*% beta[[i]][, 1:(min(nclass - 1, nz))])$class
              }
          }
      }
    }
    pred
}
