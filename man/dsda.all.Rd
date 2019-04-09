\name{dsda.all}
\alias{dsda.all}
\title{Direct sparse discriminant analysis}
\description{
Performs direct sparse discriminant analysis, with the optimal lambda chosen by cross validation. The function can perform prediction on test data as well.
}
\usage{
dsda.all(x, y, x.test.matrix=NULL, y.test=NULL, standardize=FALSE, 
lambda.opt="min", nfolds=10, lambda=lambda, alpha=1, eps=1e-7)
}
\arguments{
\item{x}{An n by p matrix containing the predictors.}
\item{y}{An n-dimensional vector containing the class labels 1 and 2.}
\item{x.test.matrix}{The predictors of a testing set. (Optional.)}
\item{y.test}{The class labels of the testing set. (Required if x.test.matrix is supplied, but otherwise optional.)}
\item{standardize}{A logic object indicating whether x.matrix should be standardized before performing DSDA. Default is FALSE.}
\item{lambda.opt}{Should be either "min" or "max", specifying whether the smallest or the largest lambda with the smallest cross validation error should be used for the final classification rule.}
\item{nfolds}{The number of folds to be used in cross validation. Default is 10.}
\item{lambda}{A sequence of lambda's.}
\item{alpha}{The elasticnet mixing parameter, the same as in glmnet. Default is alpha=1 so that the lasso penalty is used.}
\item{eps}{Convergence threshold for coordinate descent, the same as in glmnet. Default is 1e-7.}
}

\references{Mai, Q., Zou, H. and Yuan, M., (2012), "A direct approach to sparse discriminant analysis in ultra-high dimensions." Biometrika, 99, 29-42.
}

\author{
Yuqing Pan, Qing Mai, Xin Zhang
}

\seealso{
\code{\link{dsda}}
}
\value{
\item{error}{Testing error if x.test.matrix is supplied.}
\item{beta}{The coefficients of the classification rule corresponding to the optimal lambda chosen by cross validation.}
\item{s}{The optimal lambda chosen by cross validation.}
}
\examples{
data(GDS1615)   ##load the prostate data
x<-GDS1615$x
y<-GDS1615$y

x=x[which(y<3),]
y=y[which(y<3)]

n<-length(y)    ##split the original dataset to a training set and a testing set
n.test<-round(n/3)
set.seed(20120822)
id<-sample(n,n.test,replace=FALSE)
x.train<-x[-id,]
x.test<-x[id,]
y.train<-y[-id]
y.test<-y[id]

set.seed(123)
##perform direct sparse discriminant analysis
obj<-dsda.all(x.train,y.train,x.test,y.test)  
obj$error 
}