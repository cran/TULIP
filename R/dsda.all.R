dsda.all <- function(x,y,x.test.matrix=NULL,y.test=NULL,standardize=FALSE,lambda.opt="min",nfolds=10,lambda=lambda,alpha=1,eps=1e-7){
 n <- nrow(x)
 d <- ncol(x)
 n2 <- sum(y)
 n1 <- n-n2
 
 la <- cv.dsda(x, y=y, standardize=standardize, alpha=alpha, nfolds=nfolds, lambda=lambda, eps=eps, lambda.opt=lambda.opt) 
 s <- la$lambda.min
 
 obj.path <- dsda(x, y=y, lambda=s, standardize=standardize, eps=eps)
 beta <- obj.path$beta
 
 if(!missing(x.test.matrix)){
   n.test <- dim(x.test.matrix)[1]
   pred <- predict.dsda(obj.path, x.test.matrix)
   error <- mean(pred!=y.test)
   }else{error<-NA}
 result <- list(error=error, beta=beta, s=s)}
