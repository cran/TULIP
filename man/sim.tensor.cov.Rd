\name{sim.tensor.cov}
\alias{sim.tensor.cov}
\title{Simulate data}
\description{
  Simulate a data set with tensor predictor and covariates.}
\usage{
  sim.tensor.cov(tesize = 100)
}
\arguments{
\item{tesize}{Number of observations in testing data.}
}

\details{
The function simulates a data set with \eqn{10\times 10\times 10} tensor and covariate being a two-dimensional vector. Response are binary.
}
\author{
  Yuqing Pan, Qing Mai, Xin Zhang
}
\value{
  \item{x}{Simulated tensor predictor. }
  \item{z}{Simulated covariate.}
  \item{testx}{Simulated testing tensor predictor.}
  \item{testz}{Simualted testing covariate.}
  \item{vec_x}{Vectorization of \code{x}.}
  \item{vec_testx}{Vectorization of \code{testx}.}
  \item{y}{Response corresponding to \code{x} and \code{z}.}
  \item{testy}{Response corresponding to \code{testx} and \code{testz}.}
}
    