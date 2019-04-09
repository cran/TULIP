dsda_noadj<-function(x, y, standardize = FALSE, lambda = lambda, alpha = 1, eps = 1e-07){
     n<-length(y)
     n1<-sum(y==1)
     n2<-sum(y==2)
 
     obj<-glmnet(x,y,standardize=standardize,family="gaussian",alpha=alpha,thresh=eps)
     if(missing(lambda) || is.null(lambda))lambda<-obj$lambda
     beta<-coef(obj,s=lambda,standardize=standardize)
      beta1<-as.matrix(beta[-1,,drop=FALSE])
      sel<-apply(as.matrix(abs(beta1)),1,sum)
      sel<-which(sel!=0)
      mu1<-as.matrix(apply(x[y==1,sel,drop=F],2,mean))
      mu2<-as.matrix(apply(x[y==2,sel,drop=F],2,mean))
      sigma.hat<-(cov(x[y==1,sel,drop=F])*(n1-1)+cov(x[y==2,sel,drop=F])*(n2-1))/(n-2)
     beta1[sel,]<-sweep(as.matrix(beta1[sel,,drop=F]),2,sign(t(mu2-mu1)%*%beta1[sel,,drop=F]),"*")
     beta0<--t((mu1+mu2))%*%beta1[sel,,drop=F]/2-log(n1/n2)*diag(t(beta1[sel,,drop=F])%*%sigma.hat%*%beta1[sel,,drop=F])/(t(mu2-mu1)%*%beta1[sel,,drop=F])
     beta0[is.na(beta0)]<-n2-n1
     beta[1,]<-beta0

     outlist <- list(beta=beta,lambda=lambda,x=x,y=y)
     class(outlist) <- c("dsda")
     outlist
}
