predict.dsda<-function(object,newx, z=NULL, ztest=NULL, gamma=NULL, ...){

  if (is.null(z)){
    beta<-object$beta
    pred<-as.matrix(cbind(1,newx)%*%beta)
    pred<-ifelse(pred>0,2,1)
  }else{
    
    beta <- object$beta
    pred <- as.matrix(cbind(1,newx)%*%beta+ztest%*%gamma)
    print('a')
    print(pred)
    pred <- ifelse(pred>0,2,1)
  }
  
  pred

  
}
