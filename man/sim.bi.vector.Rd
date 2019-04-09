\name{sim.bi.vector}
\alias{sim.bi.vector}
\title{Simulate data}
\description{
  Simulate a binary data set with vector predictor.}
\usage{
  sim.bi.vector(tesize = 100)
}
\arguments{
\item{tesize}{Number of observations in testing data.}
}

\details{
The function simulates a data set with \eqn{p=500}. Response are binary.
}
\author{
  Yuqing Pan, Qing Mai, Xin Zhang
}
\value{
  \item{x}{Simulated vector predictor. }
  \item{testx}{Simulated testing vector predictor.}
  \item{y}{Response corresponding to \code{x}.}
  \item{testy}{Response corresponding to \code{testx}.}
}
    