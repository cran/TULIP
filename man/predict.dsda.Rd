\name{predict.dsda}
\alias{predict.dsda}
\title{Prediction for direct sparse discriminant analysis}
\description{
Predict the class labels by direct sparse discriminant analysis.
}
\usage{
\method{predict}{dsda}(object, newx, z=NULL, ztest=NULL, gamma=NULL,...)
}
\arguments{
\item{object}{An object returned by \code{dsda} or \code{msda} with binary setting.}
\item{newx}{An n by p matrix containing the predictors.}
\item{z}{Input training covariates matrix. \code{z} can be omitted if there is no covariate. }
\item{ztest}{Input testing covariates matrix. \code{ztest} can be omitted if there is no covariate.}
\item{gamma}{Coefficients of covariates obtained from \code{\link{adjvec}}. \code{gamma} is \code{NULL} if there is no covariate.}
\item{\dots}{Other arguments that can be passed to \code{predict}.}
}

\references{
Mai, Q., Zou, H. and Yuan, M. (2013), "A direct approach to sparse discriminant analysis in ultra-high dimensions." Biometrika, 99, 29-42.
}
\seealso{
\code{\link{dsda}}, \code{\link{dsda.all}}, \code{\link{predict.msda}}
}
\value{
\item{pred}{The the predicted class labels.}
}
