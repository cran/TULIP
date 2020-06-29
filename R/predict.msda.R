predict.msda <- function(object, newx, z=NULL, ztest=NULL, gamma=NULL, ...){
  switch(class(object),
         "dsda" = predict.dsda(object,newx, z=z, ztest=ztest, gamma=gamma,...),
         "msda.original" = predict.msda.original(object,newx, z=z, ztest=ztest, gamma=gamma, ...),
         "msda.modified" = predict.msda.modified(object,newx, z=z, ztest=ztest, gamma=gamma, ...),
         stop('Fitted model type is not available for prediction')
  )
}
