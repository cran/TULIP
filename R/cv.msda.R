cv.msda <- function(x, y,model=NULL,nfolds=5, lambda=NULL, lambda.opt="min",...) {
  nvars=dim(x)[2]
  nclass=length(unique(y))
  if (is.null(model)){
    if (nclass==2){
      cv.dsda(x,y,nfolds=5, lambda=NULL, lambda.opt="min",...)
    }else if (nclass>=3){
      if (nvars<=2000){ 
        cv.Omsda(x,y,nfolds=5, lambda=NULL, lambda.opt="min",...)
      }else{
        cv.Mmsda(x,y,nfolds=5, lambda=NULL, lambda.opt="min",...)
      }
    }else if (nclass==1){
      stop("y should be binary or multi-class variables.")
    }
  }else{
    switch(model,
           "binary" = cv.dsda(x,y,nfolds=5, lambda=NULL, lambda.opt="min",...),
           "multi.original" = cv.Omsda(x,y,nfolds=5, lambda=NULL, lambda.opt="min",...),
           "multi.modified" = cv.Mmsda(x,y,nfolds=5, lambda=NULL, lambda.opt="min",...),
            stop("Model should be one of 'binary', 'multi.original' and 'multi.modified'")
    )
  }
  
}