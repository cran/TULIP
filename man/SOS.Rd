\name{SOS}
\alias{SOS}
\title{Solution path for sparse discriminant analysis}
\description{
  Compute the solution path for sparse optimal scoring (SOS).}
\usage{
  SOS(x,y,standardize=FALSE,lambda=NULL,eps=1e-7)
}
\arguments{
  \item{x}{Input matrix of predictors. \code{x} is of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{y}{An n-dimensional vector containing the class labels. The classes have to be labeled as 1 and 2.}
  \item{standardize}{A logic object indicating whether x should be standardized before performing SOS. Default is FALSE.}
  \item{lambda}{A sequence of lambda's. If lambda is missed, the function will automatically generates a sequence of lambda's to fit model.}
  \item{eps}{Convergence threshold for coordinate descent, the same as in glmnet. Default is 1e-7.}
}
\references{
  Mai, Q. and Zou, H. (2013), "A note on the connection and equivalence of three sparse linear discriminant analysis methods." Technometrics, 55, 243-246.
}
\details{
The function obtains the solution path of sparse optimal scoring model through \code{\link{dsda}}.
}
\author{
  Yuqing Pan, Qing Mai, Xin Zhang
}
\value{
  \item{beta}{Output variable coefficients for each lambda. }
  \item{lambda}{The sequence of lambda's used in computing the solution path.}
  }
    \examples{
    data(GDS1615)   ##load the prostate data
    x<-GDS1615$x
    y<-GDS1615$y
    x=x[which(y<3),]
    y=y[which(y<3)]
    obj.path <- SOS(x, y)
    
    }
    