msda <- function(x, z=NULL, y, testx=NULL,testz=NULL,model=NULL, lambda = NULL, standardize = FALSE, alpha = 1, nlambda = 100,
                     lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2, 0.001), dfmax = nobs, 
                     pmax = min(dfmax * 2 + 20,  nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+06, 
                     sml = 1e-06, verbose = FALSE, perturb = NULL){

  pred = NULL
  nobs=length(y)
  nvars=dim(x)[2]
  nclass=length(unique(as.factor(y)))
  if (is.null(testx)){
    if (is.null(z)){
      objm <- msda_noadj(x,y,model=model,lambda, standardize, alpha,nlambda, lambda.factor, dfmax, pmax, pf, eps, maxit, sml, verbose, perturb)
    }else{
      obj <- adjvec(x,z,y)
      objm <- msda_noadj(obj$xres,y,model=model,lambda, standardize, alpha,nlambda, lambda.factor, dfmax, pmax, pf, eps, maxit, sml, verbose, perturb)
    }
  }
  else{  
    if (is.null(z) && is.null(testz)){
      objm <- msda_noadj(x,y,model=model,lambda, standardize, alpha,nlambda, lambda.factor, dfmax, pmax, pf, eps, maxit, sml, verbose, perturb)
      pred <- predict.msda(objm,testx)
    }else{
      if (is.null(z)){
        stop('Covariates for training data are missing.')
      }
      if (is.null(testz)){
        stop('Covariates for testing data are missing.')
      }
      obj <- adjvec(x,z,y,testx,testz)
      objm <- msda_noadj(obj$xres,y,model=model,lambda, standardize, alpha,nlambda, lambda.factor, dfmax, pmax, pf, eps, maxit, sml, verbose, perturb)
      pred <- predict.msda(objm,obj$testxres,z,testz,obj$gamma)
    }
  }
  objm$pred=pred
  outlist=c(objm)
  class(outlist) <- class(objm)
  return(outlist)
  
}

