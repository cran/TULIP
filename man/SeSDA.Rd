\name{SeSDA}
\alias{SeSDA}
\title{Solution path for semiparametric sparse discriminant analysis}
\description{
Compute the solution path for semiparametric sparse discriminant analysis.}
\usage{
SeSDA(x,y,standardize=FALSE,lambda=NULL,alpha=1,eps=1e-7)
}
\arguments{
\item{x}{Input matrix of predictors. \code{x} is of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
\item{y}{An n-dimensional vector containing the class labels. The classes have to be labeled as 1 and 2.}
\item{standardize}{A logic object indicating whether x should be standardized after transformation but before fitting classifier. Default is FALSE.}
\item{lambda}{A sequence of lambda's. If lambda is missed or NULL, the function will automatically generates a sequence of lambda's to fit model.}
\item{alpha}{The elasticnet mixing parameter, the same as in glmnet. Default is alpha=1 so that the lasso penalty is used.}
\item{eps}{Convergence threshold for coordinate descent, the same as in glmnet. Default is 1e-7.}
}
\references{
Mai, Q., Zou, H. and Yuan, M. (2013). A direct approach to sparse discriminant analysis in ultra-high dimensions. Biometrika, 99, 29-42.
}
\author{
Yuqing Pan, Qing Mai, Xin Zhang
}
\value{
\item{transform}{The tranformation functions.}
\item{objdsda}{A DSDA object fitted on transformed data.}
}
\examples{
  data(GDS1615)   ##load the prostate data
  x<-GDS1615$x
  y<-GDS1615$y
  x=x[which(y<3),]
  y=y[which(y<3)]
  obj.path <- SeSDA(x,y)

}
