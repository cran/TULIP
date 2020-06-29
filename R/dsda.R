dsda <- function(x, z=NULL, y, testx=NULL,testz=NULL, standardize = FALSE, lambda = lambda, alpha = 1, eps = 1e-07){
  pred <- NULL
  if (is.null(testx)){
    if (is.null(z)){
      objm <- dsda_noadj(x, y, standardize = standardize, lambda = lambda, alpha = alpha, eps = eps)
    }else{
      obj <- adjvec(x,z,y)
      objm <- dsda_noadj(obj$xres,y, standardize = standardize, lambda = lambda, alpha = alpha, eps = eps)
    }
  }
  else{  
    if (is.null(z) && is.null(testz)){
      objm <- dsda_noadj(x, y, standardize = standardize, lambda = lambda, alpha = alpha, eps = eps)
      pred <- predict.dsda(objm,testx)
    }else{
      if (is.null(z)){
        stop('Covariates for training data are missing.')
      }
      if (is.null(testz)){
        stop('Covariates for testing data are missing.')
      }
      obj <- adjvec(x,z,y,testx,testz)
      objm <- dsda_noadj(obj$xres, y, standardize = standardize, lambda = lambda, alpha = alpha, eps = eps)
      pred <- predict.dsda(objm,obj$testxres,z,testz,obj$gamma)
    }
  }
  outlist <- c(objm,list(pred=pred))
  class(outlist) <- c('dsda')
  return(outlist)
}
