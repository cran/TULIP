\name{cv.SeSDA}
\alias{cv.SeSDA}
\title{Cross validation for semiparametric sparse discriminant analysis}
\description{
Choose the optimal lambda for semiparametric sparse discriminant analysis by cross validation.
}
\usage{
cv.SeSDA(x, y, nfolds = 5, lambda=NULL, lambda.opt="min",
  standardize=FALSE, alpha=1, eps=1e-7)
}
\arguments{
\item{x}{An n by p matrix containing the predictors.}
\item{y}{An n-dimensional vector containing the class labels.}
\item{nfolds}{The number of folds to be used in cross validation. Default is 5.}
\item{lambda}{A sequence of lambda's.}
\item{lambda.opt}{Should be either "min" or "max", specifying whether the smallest or the largest lambda with the smallest cross validation error should be used for the final classification rule.}
\item{standardize}{A logic object indicating whether x.matrix should be standardized before performing DSDA. Default is FALSE.}
\item{alpha}{The elasticnet mixing parameter, the same as in glmnet. Default is alpha=1 so that the lasso penalty is used.}
\item{eps}{Convergence threshold for coordinate descent, the same as in glmnet. Default is 1e-7.}
}

\references{
Mai, Q., Zou, H. and Yuan, M. (2013). A direct approach to sparse discriminant analysis in ultra-high dimensions. Biometrika, 99, 29-42.
}

\seealso{
\code{cv.dsda}
\code{SeSDA}
}
\value{
\item{transform}{The transformation functions.}
\item{objdsda}{The output of cross validation from \code{cv.dsda} on transformed data.}
}
