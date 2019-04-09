\name{adjvec}
\alias{adjvec}
\title{Adjust vector for covariates.
}
\description{
Adjusts vector with respect to covariates. Vector depends on the covariates through a linear regression model. The function returns the coefficients of covariates in regression and adjusted predictor matrix for further classifier modeling. It estimates coefficients based on training data, and then adjusts training tensor. When testing data is provided, the function will automatically adjust testing data by learned coefficients as well.
}
\usage{
adjvec(x, z, y, testx = NULL, testz = NULL, is.centered = FALSE)
}
\arguments{
  \item{x}{Input matrix of dimension \eqn{N\times p}{N*p}, where \eqn{N}{N} is the number of observations and \eqn{p} is the number of variables. Each row is an observation}
  \item{z}{Input covariate matrix of dimension \eqn{N \times q}{N*q}, where \eqn{q<N}{q<N}. Each row of \code{z} is an observation.}
  \item{y}{Class label vector of dimention \eqn{N\times 1}{N*1}. For \code{K} class problems, \code{y} takes values in \eqn{\{1,\cdots,\code{K}\}}.}
  \item{testx}{Input testing matrix. Each row is a test case. When \code{testx} is not provided, the function will only adjust training data.}
  \item{testz}{Input testing covariate matrix with each row being an observation.}
  \item{is.centered}{Indicates whether the input vector and covariates have already been centered by their within class mean or not. If \code{is.centered} is FALSE, the function \code{adjvec} will center data by class. If \code{is.centered} is TRUE, the function will skip the centering step.}
}
\details{
Similar as CATCH model, assume the linear relationship between vector predictors and covariates as
\deqn{\mathbf{X}=\boldsymbol{\mu}_k+\boldsymbol{\alpha}\times\mathbf{Z}+\mathbf{E},} where \eqn{\mathbf{X}} is a \eqn{N\times p} matrix and \eqn{\boldsymbol{\alpha}} is the matrix of estimated coefficient of covariates.
The function removes the effects of covariates on response variable through vector and obtain \eqn{\mathbf{X}-\boldsymbol{\alpha}\times\mathbf{Z}} as adjusted predictors to fit MSDA and DSDA model.
}
\value{
  \item{gamma}{The estimated coefficients of covariates to plug in classifier. \code{gamma} is similar as the \eqn{\boldsymbol{\gamma}_k} defined function \code{\link{catch}} of dimension \eqn{q\times (K-1)}{q*(K-1)}, where \code{q} is the size of covariates and \code{K} is the number of classes.}
  \item{xres}{Adjusted training predictor matrix \eqn{\mathbf{X}-\boldsymbol{\alpha}\times\mathbf{Z}} after adjusting for covariates. The effect of the covariate is removed.}
  \item{testxres}{Adjusted testing predictor matrix \eqn{\mathbf{X}-\boldsymbol{\alpha}\times\mathbf{Z}} after adjusting for covariates. The effect of the covariate is removed.}
}
\references{
Pan, Y., Mai, Q., and Zhang, X. (2018), "Covariate-Adjusted Tensor Classification in High-Dimensions." Journal of the American Statistical Association, \emph{accepted}.
}
\author{
Yuqing Pan, Qing Mai, Xin Zhang
}

\seealso{
\code{\link{adjten}}
}
\examples{
n <- 50
p <- 200
k <- 2
q <- 2
x <- matrix(rnorm(n*p), n, p)
z <- matrix(rnorm(n*q), n, q)
x[1:20, ] <- x[1:20, ] + 2
z[1:20, ] <- z[1:20, ] + 0.5
y <- c(rep(1, 20), rep(2, 30))
obj <- adjvec(x, z, y)
}


