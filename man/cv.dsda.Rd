\name{cv.dsda}
\alias{cv.dsda}
\title{Cross validation for direct sparse discriminant analysis}
\description{
Choose the optimal lambda for direct sparse discriminant analysis by cross validation.
}
\usage{
cv.dsda(x, y, nfolds = 5, lambda=lambda, lambda.opt="min", 
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
\code{predict.dsda}
\code{dsda}
}
\value{
\item{lambda}{The sequence of lambda's used in cross validation.}
\item{cvm}{Cross validation errors.}
\item{cvsd}{The standard error of the cross validation errors.}
\item{lambda.min}{The optimal lambda chosen by cross validation.}
\item{model.fit}{The fitted model.}
}
