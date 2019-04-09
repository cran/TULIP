\name{predict.SeSDA}
\alias{predict.SeSDA}
\title{Prediction for semiparametric sparse discriminant analysis}
\description{
Predict the class labels by semiparametric sparse discriminant analysis.
}
\usage{
\method{predict}{SeSDA}(object, x.test,...)
}
\arguments{
\item{object}{An object returned by \code{SeSDA}.}
\item{x.test}{An n by p matrix containing the predictors.}
\item{\dots}{Other arguments that can be passed to \code{predict}.}
}

\references{
Mai, Q., Zou, H. and Yuan, M. (2013), "A direct approach to sparse discriminant analysis in ultra-high dimensions." Biometrika, 99, 29-42.
}
\seealso{
\code{\link{dsda}}, \code{\link{SeSDA}}
}
\value{
\item{pred}{The the predicted class labels.}
}
