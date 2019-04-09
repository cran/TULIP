predict.SeSDA <- function(object, x.test,...){
  n <- nrow(x.test)
  p <- ncol(x.test)
  x.test.norm <- matrix(0,n,p)
  for (i in 1:p){
    x.test.norm[,i] <- as.matrix(object$transform[[i]](x.test[,i]))
  }
  pred <- predict.dsda(object$objdsda, x.test.norm,...)
  pred

}