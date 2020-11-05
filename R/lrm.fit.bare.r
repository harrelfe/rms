##' Bare Bones Logistic Regression Fit
##'
##' This is a stripped down version of the `lrm.fit()` function that computes only the regression coefficients, variance-covariance-matrix, and log likelihood (for null and fitted model) and does not compute any model fit indexes etc.  This is for speed in simulations or with bootstrapping.  Missing data are not allowed.  The function handles binary and ordinal logistic regression (proportional odds model).
##' @title lrm.fit.bare
##' @param x a vector of matrix of covariate values
##' @param y a numeric or factor vector representing the dependent variable
##' @param maxit maximum number of iteractions
##' @param eps stopping criterion (change in -2 log likelihood)
##' @param tol matrix inversion tolerance for singularities
##' @return a list with elements `coefficients`, `var`, `fail`, `freq`, `deviance`
##' @author Frank Harrell
##' @md
lrm.fit.bare <- function(x, y, maxit=12, eps=.025, tol=1E-7) {
  opts <- double(12)
  opts[1:3] <- c(tol, eps, maxit)

  n <- length(y)
  if(n < 3) stop("must have >=3 non-missing observations")
  if(! is.matrix(x)) x <- as.matrix(x)
  dx <- dim(x)
  nx <- dx[2]
  if(nx == 0) stop('must have at least one x')
  if(dx[1] != n) stop("x and y must have same length")
  storage.mode(x) <- "double"

  est <- 1 : nx
  xname <- dimnames(x)[[2]]
  if(! length(xname)) xname <- paste("x[", 1 : nx, "]", sep="")

  if(! is.factor(y)) y <- as.factor(y)
  y       <- unclass(y)
  ylevels <- levels(y)

  kint   <- as.integer(length(ylevels) - 1)
  ftable <- integer(5001 * (kint + 1))
  numy        <- tabulate(y)
  names(numy) <- ylevels
  y   <- as.integer(y - 1)
  nvi <- as.integer(nx + kint)

  weights <- rep(1., n)
  storage.mode(weights) <- 'double'
  
  sumwty <- tapply(weights, y, sum)
  sumwt  <- sum(sumwty)
  sumw <- as.integer(round(sumwty))
  
  ncum <- rev(cumsum(rev(sumwty)))[2 : (kint + 1)]
  pp   <- ncum / sumwt
  initial <- rep(0., nvi)
  initial[1 : kint] <- log(pp / (1 - pp))
  storage.mode(initial) <- "double"
  
  loglik <- -2 * sum(sumwty * logb(sumwty / sum(sumwty)))
  ## loglik <-  -2 * sum(numy * logb(numy/n))
  
  penmat <- matrix(0, ncol=nvi, nrow=nvi)
  storage.mode(penmat) <- 'double'

  z <- 
    .Fortran(F_lrmfit, coef=initial, nx, est, x, y, offset=0.,
             u=double(nvi),
             double(nvi * (nvi + 1) / 2), loglik=double(1), n, nx, sumw, nvi,
             v=double(nvi * nvi), double(nvi), double(2 * nvi), double(nvi),
             pivot=integer(nvi), opts=opts, ftable, penmat, weights)

	irank <- z$opts[7]
	if(irank < nvi) {
      cat("singular information matrix in lrm.fit.bare (rank=", irank,
          ").  Offending variable(s):\n")
      cat(paste(xname[est[z$pivot[nvi : (irank + 1)] - kint]],
                collapse=" "), "\n")
      return(structure(list(fail=TRUE)))
    }
	loglik <- c(loglik, z$loglik)
  
  dvrg <- z$opts[6] > 0
  
  
  ## Invert v with respect to fitted variables
  info.matrix <- matrix(z$v, nrow=nvi, ncol=nvi)
  v <- solvet(info.matrix, tol=tol)
  irank <- nvi

  name <- if(kint == 1) "Intercept"
          else 
            paste("y>=", ylevels[2 : (kint + 1)], sep="")
  name <- c(name, xname)
  kof <- z$coef

  names(kof) <- name
  dimnames(v) <- list(name, name)
 
  list(freq=numy, fail=dvrg, coefficients=kof,
       var=v, deviance=loglik)
}
