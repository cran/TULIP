\name{cv.msda}
\alias{cv.msda}
\title{
Cross-validation for DSDA/MSDA through function \code{msda}
}
\description{
Performs K-fold cross validation for \code{msda} and returns the best tuning parameter \eqn{\lambda} in the user-specified or automatically generated choices.
}
\usage{
cv.msda(x, y, model = NULL, nfolds = 5, lambda = NULL,
 lambda.opt = "min", ...)
}
\arguments{
  \item{x}{Input matrix of predictors. \code{x} is of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{y}{Class label. For \code{K} class problems, \code{y} takes values in \eqn{\{1,\cdots,\code{K}\}}.}
  \item{model}{Method type. The \code{model} argument can be one of \code{'binary'}, \code{'multi.original'}, \code{'multi.modified'} and the default is NULL. The function supports fitting DSDA and MSDA models by specifying method type. Without specification, the function will automatically choose one of the methods. If the response variable is binary, the function will fit a DSDA model. If the response variable is multi-class, the function will fit an original MSDA model for dimension \eqn{p<=2000} and a modified MSDA model for dimension \eqn{p>2000}.}
  \item{nfolds}{
	Number of folds. Default value is 5. Although \code{nfolds}
	can be as large as the sample size (leave-one-out CV), it is not
	recommended for large datasets. Smallest value allowable is \code{nfolds=3} for \code{multi.original} and \code{multi.modified}.}
  \item{lambda}{User-specified \code{lambda} sequence for cross validation. If not specified, the algorithm will generate a sequence of \code{lambda}s based on all data and cross validate on the sequence.}
  \item{lambda.opt}{The optimal criteria when multiple elements in \code{lambda} return the same minimum classification error. "\code{min}" will return the smallest \code{lambda} with minimum cross validation error. "\code{max}" will return the largest \code{lambda} with the minimum cross validation error.
}
  \item{\dots}{
other arguments that can be passed to \code{msda}.}
}
\details{
The function \code{\link{cv.msda}} runs function \code{\link{msda}} \code{nfolds+1} times. The first one fits model on all data. If \code{lambda} is specified, it will check if all \code{lambda} satisfies the constraints of \code{dfmax} and \code{pmax} in \code{\link{msda}}. If not, a \code{lambda} sequence will be generated according to \code{lambda.factor} in \code{\link{msda}}. Then the rest \code{nfolds} many replicates will fit model on \code{nfolds-1} many folds data and predict on the omitted fold, repectively. Return the \code{lambda} with minimum average cross validation error and the largest \code{lambda} within one standard error of the minimum.

Similar as \code{\link{msda}}, user can specify which method to use by inputing argument \code{model}. Without specification, the function can automatically decide the method by number of classes and variables.
}
\value{
An object of class \code{cv.dsda} or \code{cv.msda.original} or \code{cv.msda.modified} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{The actual \code{lambda} sequence used. The user specified sequence or automatically generated sequence could be truncated by constraints on \code{dfmax} and \code{pmax}.}
  \item{cvm}{The mean of cross validation errors for each \code{lambda}.}
  \item{cvsd}{The standard error of cross validaiton errors for each \code{lambda}.}
		\item{lambda.min}{The \code{lambda} with minimum cross validation error. If \code{lambda.opt} is \code{min}, then returns the smallest \code{lambda} with minimum cross validation error. If \code{lambda.opt} is \code{max}, then returns the largest \code{lambda} with minimum cross validation error.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is
		within one standard error of the minimum. This arguement is only available for object \code{cv.msda.original} and \code{cv.msda.modified}.}
		\item{model.fit}{A fitted \code{cv.dsda} or \code{cv.msda.original} or \code{cv.msda.modified} object for the full data.}
}
\references{
Mai, Q., Zou, H. and Yuan, M. (2012), "A direct approach to sparse discriminant analysis in ultra-high dimensions." Biometrica, 99, 29-42.

Mai, Q., Yang, Y., and Zou, H. (2017), "Multiclass sparse discriminant analysis." Statistica Sinica, in press.

URL: \url{https://github.com/emeryyi/msda}\cr
}
\author{
Yuqing Pan, Qing Mai, Xin Zhang}

\seealso{
\code{\link{msda}}}
\examples{
data(GDS1615)
x <- GDS1615$x
y <- GDS1615$y
obj.cv <- cv.msda(x=x, y=y, nfolds=5, lambda.opt="max")
lambda.min <- obj.cv$lambda.min
obj <- msda(x=x, y=y, lambda=lambda.min)
pred <- predict(obj,x)
}

