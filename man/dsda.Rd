\name{dsda}
\alias{dsda}
\title{Solution path for direct sparse discriminant analysis}
\description{
Compute the solution path for direct sparse discriminant analysis (DSDA).}
\usage{
dsda(x, z=NULL, y, testx=NULL, testz=NULL, standardize=FALSE, 
 lambda=lambda, alpha=1, eps=1e-7)
}
\arguments{
\item{x}{Input matrix of predictors. \code{x} is of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
 \item{z}{Input covariate matrix of dimension \eqn{N \times q}{N*q}, where \eqn{q<N}{q<N}. \code{z} can be omitted if covariate is absent. }
\item{y}{An n-dimensional vector containing the class labels. The classes have to be labeled as 1 and 2.}

 \item{testx}{Input testing matrix. Each row is a test case. When \code{testx} is not provided, the function will only fit the model and return the classifier. When \code{testx} is provided, the function will predict response on \code{testx} as well.}
  \item{testz}{Input testing covariate matrix. Can be omitted if covariate is absent. However, training covariates \code{z} and testing covariates \code{testz} must be provided or not at the same time.}

\item{standardize}{A logic object indicating whether x should be standardized before performing DSDA. Default is FALSE.}

\item{lambda}{A sequence of lambda's. If lambda is missed, the function will automatically generates a sequence of lambda's to fit model.}
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
\item{beta}{Output variable coefficients for each lambda. The first element of each solution is the intercept.}
\item{lambda}{The sequence of lambda's used in computing the solution path.}
\item{x}{The predictor matrix in training data.}
\item{y}{The class label in training data.}
\item{pred}{Predicted categorical response for each value in sequence \code{lambda} when \code{testx} is provided.}
}
\examples{
  data(GDS1615)   ##load the prostate data
  x<-GDS1615$x
  y<-GDS1615$y
  x=x[which(y<3),]
  y=y[which(y<3)]
  obj.path <- dsda(x, y=y)

}
