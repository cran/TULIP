 cv.dsda <- function (x, y, nfolds = 5,lambda=lambda, lambda.opt="min", standardize=FALSE, alpha=1, eps=1e-7)
 {
     n <- length(y)
     n1 <- sum(y==1)
     n2 <- sum(y==2)

     if (nfolds>n) stop("The number of folds should be smaller than the sample size.")

     all.folds <- cv.folds(length(y), nfolds)
     
     if(missing(lambda)||is.null(lambda)){
        fit <- glmnet(x, y,  family="gaussian",alpha=alpha,standardize=standardize,thresh=eps)
        lambda <- fit$lambda}

     nlambda <- length(lambda)
     residmat <- matrix(0, nlambda, nfolds)
     for (i in seq(nfolds)) {
         omit <- all.folds[[i]]
         fit <- dsda_noadj(x[-omit, , drop = FALSE], y[-omit], lambda=lambda, standardize=standardize, alpha=alpha, eps=eps)
         fit <- predict.dsda(fit,x[omit,,drop=FALSE])
         residmat[, i] <- apply(abs(sweep(fit,1,y[omit])),2,mean)}
       
     residmat[is.na(residmat)] <- min(n1/n,n2/n)
     residmat <- matrix(residmat,nrow=nlambda)
     cv <- apply(residmat, 1, mean)
     cv.error <- sqrt(apply(residmat, 1, var)/nfolds)
     if(lambda.opt=="min"){
        bestlambda <- min(lambda[which(cv==min(cv))])}
     else{
        bestlambda <- max(lambda[which(cv==min(cv))])}
     object <- list(lambda = lambda, cvm = cv, cvsd = cv.error,lambda.min=bestlambda, model.fit=fit)
     class(object) <- c("cv.dsda")
     invisible(object)
 }
 
