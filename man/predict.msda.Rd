\name{predict.msda}
\alias{predict.msda}
\title{
Predict categorical responses for vector data.
}
\description{
Predict categorical responses on new vector data given the fitted DSDA/MSDA model input.
}
\usage{
\method{predict}{msda}(object, newx, z = NULL, ztest = NULL, gamma = NULL,...)
}

\arguments{
\item{object}{Fitted model object from \code{msda}. The model object can be anyone of \code{binary}, \code{multi.original} and \code{multi.modified}.}
\item{newx}{The matrix of new values for \code{x} at which predictions are
to be made. If covariates exist, then \code{newx} should be adjusted matrix.}
\item{z}{Input training covariates matrix. \code{z} can be omitted if there is no covariate. }
\item{ztest}{Input testing covariates matrix. \code{ztest} can be omitted if there is no covariate.}
\item{gamma}{Coefficients of covariates obtained from \code{\link{adjvec}}. \code{gamma} is \code{NULL} if there is no covariate.}

\item{\dots}{Other arguments that can be passed to \code{predict}. }
}

\value{
Predicted class label(s) at the entire sequence of the penalty parameter \code{lambda} used to create the model.
}
\details{
The function fits LDA model on selected discriminant vectors. Call \code{\link{predict}} or \code{\link{predict.msda}} to perform prediction. When covariates exist, users could first call \code{\link{adjvec}} to make adjustment and obtain obtain \code{gamma}. The fitted model from \code{\link{msda}} should also takes adjusted vector as input. The \code{newx} in \code{\link{predict.msda}} shoudl be adjusted vector as well.
}

\references{
Mai, Q., Zou, H. and Yuan, M. (2012), "A direct approach to sparse discriminant analysis in ultra-high dimensions." Biometrica, 99, 29-42.

Mai, Q., Yang, Y., and Zou, H. (2017), "Multiclass sparse discriminant analysis." Statistica Sinica, in press.

Pan, Y., Mai, Q., and Zhang, X. (2018), "Covariate-Adjusted Tensor Classification in High-Dimensions." Journal of the American Statistical Association, \emph{accepted}.

}
\author{
Yuqing Pan, Qing Mai, Xin Zhang
}

\seealso{
\code{\link{msda}}}
\examples{
data(GDS1615)
x<-GDS1615$x
y<-GDS1615$y
obj <- msda(x = x, y = y)
pred<-predict(obj,x)
}

