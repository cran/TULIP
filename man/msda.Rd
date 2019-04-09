\name{msda}
\alias{msda}
\title{
Fits a regularization path of Sparse Discriminant Analysis and predicts}
\description{
Fits a regularization path of Sparse Discriminant Analysis at a sequence of regularization parameters lambda. Performs prediction when testing data is provided. The \code{msda} function solves classification problem by fitting a sparse discriminant analysis model. When covariates are provided, the function will first make adjustment on the training data. It provides three models: \code{binary} for fitting DSDA model to solve binary classification problems, \code{multi.original} and \code{multi.modified} for fitting MSDA model to solve multi-class classification problems. \code{multi.original} runs faster for small dimension case but the computation ability is limited to a relatively large dimension. \code{multi.modified} has no such limitation and works in ultra-high dimensions. User can specify method by argument or use the default settings.
}
\usage{
msda(x, z=NULL, y, testx=NULL,testz=NULL, model = NULL, lambda = NULL, 
 standardize=FALSE, alpha=1, nlambda = 100, 
 lambda.factor = ifelse((nobs - nclass)<= nvars, 0.2, 1e-03), dfmax = nobs, 
 pmax = min(dfmax * 2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, 
 maxit = 1e+06, sml = 1e-06, verbose = FALSE, perturb = NULL)
}
\arguments{
  \item{x}{Input matrix of predictors. \code{x} is of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{z}{Input covariate matrix of dimension \eqn{N \times q}{N*q}, where \eqn{q<N}{q<N}. \code{z} can be omitted if covariate is absent. }
  \item{y}{Class labl. This argument should be a factor for classification. For \code{model}=\code{'binary'}, \code{y} should be a binary variable with values 1 and 2. For \code{model}=\code{'multi.original'} or \code{'multi.modified'}, \code{y} should be a multi-class variable starting from 1.
}
  \item{testx}{Input testing matrix. Each row is a test case. When \code{testx} is not provided, the function will only fit the model and return the classifier. When \code{testx} is provided, the function will predict response on \code{testx} as well.}
  \item{testz}{Input testing covariate matrix. Can be omitted if covariate is absent. However, training covariates \code{z} and testing covariates \code{testz} must be provided or not at the same time.}
  \item{model}{Method type. The \code{model} argument can be one of \code{'binary'}, \code{'multi.original'}, \code{'multi.modified'} and the default is NULL. The function supports fitting DSDA and MSDA models by specifying method type. Without specification, the function will automatically choose one of the methods. If the response variable is binary, the function will fit a DSDA model. If the response variable is multi-class, the function will fit an original MSDA model for dimension \eqn{p<=2000} and a modified MSDA model for dimension \eqn{p>2000}.}
  \item{lambda}{
	A user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have the program compute its own \code{lambda} sequence based on
	\code{nlambda} and \code{lambda.factor}. Supplying a value of
	\code{lambda} overrides this. It is better to supply
	a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}

\item{standardize}{A logic object indicating whether x should be standardized before performing DSDA. Default is FALSE. This argument is only valid for \code{model = 'binary'}.}

\item{alpha}{The elasticnet mixing parameter, the same as in glmnet. Default is alpha=1 so that the lasso penalty is used in DSDA. This argument is only valid for \code{model = 'binary'}.}

  \item{nlambda}{The number of tuning values in sequence \code{lambda}. If users do not specify \code{lambda} values, the package will generate a solution path containing \code{nlambda} many tuning values of \code{lambda}. Default is 100 for \code{model = 'multi.original'} and 50 for \code{model = 'multi.modified'}.}
  \item{lambda.factor}{
	The factor for getting the minimal lambda in \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on \eqn{p}{p} (the number of predictors) and its relationship with \eqn{N} (the number of rows in the matrix of predictors). 
	For Original MSDA, if \eqn{N > p}, the default is \code{0.0001}, close to zero.  If \eqn{N<p}, the default is \code{0.2}. For Modified MSDA, if \eqn{p\le 5000}, the default is \code{0.2}. If \eqn{5000<p\le 30000}, the default is \code{0.4}. If \eqn{p>30000}, the default is \code{0.5}.
	A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is user-defined \code{lambda} sequence. This argument is only valid for \code{multi.original} and \code{multi.modified}.}

  \item{dfmax}{
	The maximum number of selected variables in the model. Default is the number of observations \code{N}. This argument is only valid for \code{multi.original} and \code{multi.modified}.}
  \item{pmax}{
	The maximum number of potential selected variables during iteration. In middle step, the algorithm can select at most \code{pmax} variables and then shrink part of them such that the nubmer of final selected variables is less than \code{dfmax}. Default is \eqn{\min(dfmax\times 2+20, N)}{min (dfmax*2+20, p)}.}
  \item{pf}{
	L1 penalty factor of length \eqn{p}{p}. Separate L1 penalty weights can be applied to each coefficient of \eqn{\theta}{theta} to allow
	differential L1 shrinkage. Can be 0 for some variables, which implies
	no L1 shrinkage, and results in that variable always being included in the
	model. Default is 1 for all variables (and implicitly infinity for
	variables listed in \code{exclude}). This argument is only valid for \code{multi.original} and \code{multi.modified}.}
  \item{eps}{
	Convergence threshold for coordinate descent. Each inner
	coordinate descent loop continues until the relative change in any
	coefficient. Defaults value is \code{1e-4}.}
  \item{maxit}{
Maximum number of outer-loop iterations allowed at fixed lambda value. Default is 1e6. If models do not converge, consider increasing \code{maxit}. This argument is only valid for \code{multi.original} and \code{multi.modified}.}
  \item{sml}{
    Threshold for ratio of loss function change after each iteration to old loss function value. Default is \code{1e-06}. This argument is only valid for \code{multi.original} and \code{multi.modified}.
}
  \item{verbose}{
Whether to print out computation progress. The default is \code{FALSE}. This argument is only valid for \code{multi.original} and \code{multi.modified}.}
  \item{perturb}{A scalar number. If it is specified, the number will be added to each diagonal element of the covariance matrix as perturbation. The default is \code{NULL}. This argument is only valid for \code{multi.original} and \code{multi.modified}.}

}
\details{

The \code{msda} function fits a linear discriminant analysis model for vector \eqn{X} as follows:
\deqn{\mathbf{X}|Y=k\sim N(\boldsymbol{\mu}_k,\boldsymbol{\Sigma}).}
The categorical response is predicted from the Bayes rule:
\deqn{\widehat{Y}=\arg\max_{k=1,\cdots,K}{(\mathbf{X}-\frac{\boldsymbol{\mu}_k}{2})^T\boldsymbol{\beta}_k+\log\pi_k}.}
The parameter \code{model} specifies which method to use in estimating \eqn{\boldsymbol{\beta}}. Users can use \code{binary} for binary problems and \code{binary} and \code{multi.modified} for multi-class problems. In \code{multi.original}, the algorithm first computes and stores \eqn{\boldsymbol{\Sigma}}, while it doesn't compute or store the entire covariance matrix in \code{multi.modified}. Since the algorithm is element-wise based, \code{multi.modified} computes each element of covariance matrix when needed. Therefore, \code{multi.original} is faster for low dimension but \code{multi.modified} can fit model for a much higher dimension case.

Note that for computing speed reason, if models are not converging or running slow, consider increasing \code{eps} and \code{sml}, or decreasing
\code{nlambda}, or increasing \code{lambda.factor} before increasing
\code{maxit}. Users can also reduce \code{dfmax} to limit the maximum number of variables in the model.

The arguments list out all parameters in the three models, but not all of them are necessary in applying one of the methods. See the specific explaination of each argument for more detail. Meanwhile, the output of DSDA model only includes \code{beta} and \code{lambda}.
}
\value{
An object with S3 class \code{dsda} or \code{msda.original} and \code{msda.modified}.
		\item{beta}{Output variable coefficients for each \code{lambda}, which is the estimation of \eqn{\boldsymbol{\beta}} in the Bayes rule. \code{beta} is a list of length being the number of \code{lambda}s. Each element of \code{beta} is a matrix of dimension \eqn{nvars\times (nclass-1)}{nvars*(nclass-1)}. For \code{model = 'dsda'}, \code{beta} is a vector of length \eqn{nvars+1}, where the first element is intercept.}
		\item{df}{The number of nonzero coefficients for each value of \code{lambda}.}
		\item{obj}{The fitted value of the objective function for each value of \code{lambda}.}
		\item{dim}{Dimension of each coefficient matrix.}
		\item{lambda}{The actual \code{lambda} sequence used. The user specified sequence or automatically generated sequence could be truncated by constraints on \code{dfmax} and \code{pmax}.}
		\item{x}{The input matrix of predictors for training.}
		\item{y}{Class label in training data.}
		\item{npasses}{Total number of iterations (the most inner loop) summed over all lambda values}
		\item{jerr}{Error flag, for warnings and errors, 0 if no error.}
		\item{sigma}{Estimated sigma matrix. This argument is only available in object \code{msda.original}.}
		\item{delta}{Estimated delta matrix. delta[k] = mu[k]-mu[1].}
		\item{mu}{Estimated mu vector.}
		\item{prior}{Prior probability that y belong to class k, estimated by mean(y that belong to k).}
		\item{call}{The call that produced this object}
    \item{pred}{Predicted categorical response for each value in sequence \code{lambda} when \code{testx} is provided.}
		
}
\references{
Mai, Q., Zou, H. and Yuan, M. (2012), "A direct approach to sparse discriminant analysis in ultra-high dimensions." Biometrica, 99, 29-42.

Mai, Q., Yang, Y., and Zou, H. (2017), "Multiclass sparse discriminant analysis." Statistica Sinica, in press.

URL: \url{https://github.com/emeryyi/msda}\cr
}
\author{
Yuqing Pan, Qing Mai, Xin Zhang
}

\seealso{\code{\link{cv.msda}}, \code{\link{predict.msda}}
}
\examples{
data(GDS1615)
x<-GDS1615$x
y<-GDS1615$y
obj <- msda(x = x, y = y)
}

