# the main function
Mmsda <- function(x, y,lambda = NULL, standardize=NULL, alpha=NULL, nlambda = 50, lambda.factor = ifelse(nvars<=5000,0.2,ifelse(nvars>=30000,0.5,0.4)), dfmax = nobs, pmax = min(dfmax * 
    2 + 20, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+06, sml = 1e-06, 
    verbose = FALSE, perturb = 0) {
    ## data setup
    this.call <- match.call()
    tmp <- Mmsda.prep(x, y)
    delta <- as.matrix(tmp$delta)
    mu <- as.matrix(tmp$mu)
    prior <- tmp$prior
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
    nclass <- as.integer(length(unique(y)))
    vnames <- colnames(x)
    if (is.null(vnames)) 
        vnames <- paste("V", seq(nvars), sep = "")
    nk <- as.integer(dim(delta)[1])
    ## parameter setup
    if (length(pf) != nvars) 
        stop("The size of penalty factor must be same as the number of input variables")
    maxit <- as.integer(maxit)
    verbose <- as.integer(verbose)
    sml <- as.double(sml)
    pf <- as.double(pf)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax)
    pmax <- as.integer(pmax)
    ## lambda setup
    nlam <- as.integer(nlambda)
    if (is.null(lambda)) {
        if (lambda.factor >= 1) 
            stop("lambda.factor should be less than 1")
        flmin <- as.double(lambda.factor)
        ulam <- double(1)  #ulam=0 if lambda is missing
    } else {
        # flmin=1 if user define lambda
        flmin <- as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))  #lambda is declining
        nlam <- as.integer(length(lambda))
    }
    ## call Fortran core
    fit <- .Fortran("Mmsda", obj = double(nlam), nk, nvars, nobs, as.double(x), as.integer(y), as.double(t(mu)), as.double(delta),
        pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), 
        theta = double(pmax * nk * nlam), itheta = integer(pmax), ntheta = integer(nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1), perturb=as.double(perturb))
    ## output
    outlist <- formatoutput(fit, maxit, pmax, nvars, vnames, nk)
    outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr, delta = delta, mu = mu, prior = prior, call = this.call))
    if (is.null(lambda)) 
        outlist$lambda <- lamfix(outlist$lambda)
    class(outlist) <- c("msda.modified")
    outlist
}
