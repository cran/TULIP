msda_noadj <- function(x, y,model=NULL, lambda =
                 NULL, standardize = FALSE, alpha = 1, nlambda = 100,
                 lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2,
                 0.001), dfmax = nobs, pmax = min(dfmax * 2 + 20,
                 nvars), pf = rep(1, nvars), eps = 1e-04, maxit =
                 1e+06, sml = 1e-06, verbose = FALSE, perturb = NULL) {

	nvars=dim(x)[2]
	nobs=length(y)
	nclass=length(unique(y))
	if (is.null(model)){
	  if (nclass==2){
	    dsda(x,y=y,lambda=lambda, standardize=standardize, alpha=alpha)
	  }else if (nclass>=3){
  	  if (nvars<=2000){ 
	  	  Omsda(x,y=y,lambda=lambda, nlambda=nlambda, lambda.factor=lambda.factor, dfmax=dfmax, pmax=pmax, pf=pf, eps=eps, maxit=maxit, sml=sml, verbose=verbose, perturb=perturb)
	    }else{
	  	  Mmsda(x,y=y,lambda=lambda, nlambda=nlambda, lambda.factor=lambda.factor, dfmax=dfmax, pmax=pmax, pf=pf, eps=eps, maxit=maxit, sml=sml, verbose=verbose, perturb=perturb)
	    }
	  }else if (nclass==1){
	    stop("y should be binary or multi-class variables.")
	  }
	}else{
      switch(model,
             "binary" = dsda(x,y=y,lambda=lambda, standardize=standardize, alpha=alpha),
	      "multi.original" = Omsda(x,y=y,lambda=lambda, nlambda=nlambda, lambda.factor=lambda.factor, dfmax=dfmax, pmax=pmax, pf=pf, eps=eps, maxit=maxit, sml=sml, verbose=verbose, perturb=perturb),
             "multi.modified" = Mmsda(x,y=y,lambda=lambda, nlambda=nlambda, lambda.factor=lambda.factor, dfmax=dfmax, pmax=pmax, pf=pf, eps=eps, maxit=maxit, sml=sml, verbose=verbose, perturb=perturb),
             stop("Model should be one of 'binary', 'multi.original' and 'multi.modified'")
             )
	}
	
}