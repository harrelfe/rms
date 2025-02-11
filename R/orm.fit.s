orm.fit <- function(x=NULL, y,
                    family=c("logistic","probit","loglog","cloglog","cauchit"),
                    offset, initial,
                    opt_method=c('NR', 'LM'),
                    maxit=30L, eps=5e-4, gradtol=1e-3, abstol=1e10,
                    minstepsize=1e-2, tol=.Machine$double.eps, trace=FALSE,
                    penalty.matrix=NULL, weights=NULL, normwt=FALSE,
                    scale=FALSE, inclpen=TRUE, y.precision = 7,
                    compstats=TRUE, onlydata=FALSE)
{
  cal        <- match.call()
  family     <- match.arg(family)
  opt_method <- match.arg(opt_method)

  n <- NROW(y)
  if(! length(x)) {
      p     <- 0
      xname <- NULL
      x     <- 0.
    }
  else  {
      if(! is.matrix(x)) x <- as.matrix(x)
      dx <- dim(x)
      p  <- dx[2L]
      if(dx[1] != n) stop("x and y must have same number of rows")
      xname <- dimnames(x)[[2]]
      if(! length(xname)) xname <- paste("x[", 1 : p, "]", sep="")
  }

  len.penmat <- length(penalty.matrix)
  penpres    <- len.penmat && any(penalty.matrix != 0.)
  if(p == 0 && penpres) stop('may not specify penalty.matrix without predictors')
  if(penpres && any(dim(penalty.matrix) != p))
    stop(paste("penalty.matrix does not have", p, "rows and columns"))
  penmat <- if(! penpres) matrix(0e0, nrow=p, ncol=p) else penalty.matrix

  ## Extreme value type I dist = Gumbel maximum = exp(-exp(-x)) = MASS:::pgumbel
  ## Gumbel minimum = 1 - exp(-exp(x))
  families <- probabilityFamilies
  familiesDefined <- names(families)
  link  <- match(family, familiesDefined, nomatch=0)
  if(link == 0)
    stop('family must be one of ', paste(familiesDefined, collapse=' '))
  fam <- families[[family]]

  wtpres <- TRUE
  if(! length(weights)) {
    wtpres  <- FALSE
    normwt  <- FALSE
    weights <- rep(1.0, n)
  }
  if(length(weights) != n) stop('length of weights must equal length of y')
  if(normwt) weights <- weights * n / sum(weights)

  initial.there <- ! missing(initial)

  if(p > 0 && scale) {
    x      <- scale(x)
    scinfo <- attributes(x)[c('scaled:center', 'scaled:scale')]
    xbar   <- as.matrix(scinfo[[1]])
    xsd    <- as.matrix(scinfo[[2]])
    # Transform penalty matrix to new scale
    trans <- rbind(cbind(1., matrix(0., 1, p)),   # 1.: only dealing with middle intercept
                   cbind(- matrix(xbar / xsd, ncol=1),
                          diag(1. / as.vector(xsd), ncol=p)))
    if(penpres) penmat <- t(trans[-1, -1]) %*% penmat %*% trans[-1, -1]
    }


  # model.extract which calls model.response may not keep Ocens class
  isOcens <- NCOL(y) == 2
  orange  <- NULL
  yupper  <- NULL
  npsurv  <- NULL
  ctype   <- integer(n)

  if(isOcens) {
    Y       <- y                # save original values for initializing intercepts with survfit
    a       <- attributes(y)
    ylevels <- a$levels
    yupper  <- a$upper
    npsurv  <- a$npsurv
    if(missing(initial) && ! length(npsurv))
      stop('you must specify npsurv=TRUE to Ocens when not specifying initial to orm.fit')
    k       <- length(ylevels) - 1
    # Ocens uses -Inf and Inf for left and right censoring; later we transform to integer codes 0, k+1
    ctype[is.infinite(y[, 1])] <- 1               # left censoring
    ctype[is.infinite(y[, 2])] <- 2               # right
    ctype[ctype == 0 & (y[, 1] < y[, 2])] <- 3    # interval
    Ncens   <- c(left     = sum(ctype == 1),
                 right    = sum(ctype == 2),
                 interval = sum(ctype == 3) )
    y2      <- ifelse(ctype == 2, k + 1, y[, 2] - 1L)
    y       <- ifelse(ctype == 1, -1L,   y[, 1] - 1L)
    numy    <- a$freq
    mediany <- a$median
    kmid    <- which.min(abs(ylevels[-1] - mediany))
    anycens <- any(ctype > 0)
    orange  <- a$range
  } else {
  w         <- recode2integer(y, precision=y.precision)
  y         <- w$y - 1
  y2        <- y
  ylevels   <- w$ylevels
  k         <- length(ylevels) - 1
  kmid      <- max(w$whichmedian - 1L, 1L)
  numy      <- w$freq
  mediany   <- w$median
  Ncens     <- c(left=0L, right=0L, interval=0L)
  anycens   <- FALSE
  }

  intcens <- 1L * any(ctype == 3)

  if(k == 1) kmid <- 1

  ylev <- ylevels[-1L]
  if(length(yupper)) ylev <- ifelse(yupper[-1] > ylev, paste(ylev, '-', yupper[-1]), ylev)
  iname <- if(k == 1) "Intercept" else paste0("y>=", ylev)
  name  <- c(iname, xname)

  if(onlydata) return(
    if(isOcens) list(Y=cbind(y, y2), Ncens=Ncens, k=k, iname=iname)
    else list(y=y, k=k, iname=iname) )

  if(missing(offset) || ! length(offset) || (length(offset) == 1 && offset == 0.))
    offset <- rep(0., n)
  ofpres <- ! all(offset == 0.)
  if(ofpres && length(offset) != n) stop("offset and y must have same length")

  if(n < 3) stop("must have >=3 non-missing observations")

  nv <- p + k

  # These are used for initial intercept estimates when there is no censoring (OK)
  # and also in R2 statistics
  sumwty <- tapply(weights, y, sum)
  sumwt  <- sum(sumwty)
  if(anycens) {
    ess <- n - sum(Ncens)
    Nw  <- n
  } else {
    rfrq   <- sumwty / sumwt
    ess    <- n * (1 - sum(rfrq ^ 3))
    Nw     <- sumwt
  }

  if(missing(initial)) {
    if(anycens) {
      # If only censored values are right censored, then MLEs of intercepts are exactly
      # the link function of Kaplan-Meier estimates.  We're ignoring weights.
      # This should also work for only left censoring, and reasonable results
      # are expected for interval and mixed censoring using a Turnbull-type estimator
      # computed from icenReg::ic_np as stored in the npsurv object from Ocens.
      # ic_np also handles non-interval censoring.
      pp <- npsurv$surv[-1]
      # n.risk <- attr(pp, 'n.risk')
    } else {
      ncum   <- rev(cumsum(rev(sumwty)))[2 : (k + 1)]
      pp     <- ncum / sumwt
    }
    finverse <- eval(fam[2])
    initial <- finverse(pp)
    if(ofpres) initial <- initial - mean(offset)
    initial <- c(initial, rep(0., p))
 }

  if(! anycens) loglik <- -2 * sum(sumwty * log(sumwty / sum(sumwty)))

  if(anycens || (p==0 & ! ofpres)) {
    z <- ormfit(NULL, y, y2, k, intcens, initial=initial[1 : k],
                offset=offset, wt=weights,
                penmat=penmat, opt_method=opt_method, maxit=maxit,
                tolsolve=tol, objtol=eps, gradtol=gradtol, paramtol=abstol,
                trace=trace, link=link, iname=iname, xname=xname)
    if(z$fail) return(structure(list(fail=TRUE), class="orm"))
    kof    <- z$coef
    loglik <- z$loglik
    initial <- c(kof, rep(0., p))
    info   <- z$info
  }

  if(ofpres) {
    ## Fit model with only intercept(s) and offset
    ## Check that lrm.fit uses penmat in this context   ??
    z <- ormfit(NULL, y, y2, k, intcens, initial=initial[1 : k],
                offset=offset, wt=weights,
                penmat=penmat, opt_method=opt_method, maxit=maxit,
                tolsolve=tol, objtol=eps, gradtol=gradtol, paramtol=abstol,
                trace=trace, link=link, iname=iname, xname=xname)
    if(z$fail) return(structure(list(fail=TRUE), class="orm"))
    kof    <- z$coef
    loglik <- c(loglik, z$loglik)
    initial <- c(z$coef, rep(0., p))
    if(p == 0) info <- z$info
  }

  if(p > 0) {
    # Fit model with intercept(s), offset, covariables
    z <- ormfit(x, y, y2, k, intcens, initial=initial, offset=offset, wt=weights,
                penmat=penmat, opt_method=opt_method,
                maxit=maxit, tolsolve=tol, objtol=eps,
                gradtol=gradtol, paramtol=abstol,
                trace=trace, link=link, iname=iname, xname=xname)
    if(z$fail) return(structure(list(fail=TRUE), class="orm"))
    loglik <- c(loglik, z$loglik)
    kof  <- z$coef
    info <- z$info
    # Compute linear predictor before unscaling beta, as x is scaled
    lp <- matxv(x, kof, kint=kmid)

    if(scale) {
      betas         <- kof[- (1 : k)]
      kof[1 : k]    <- kof[1 : k] - sum(betas * xbar / xsd)
      kof[-(1 : k)] <- betas / xsd
      xbar          <- as.vector(xbar)
      names(xbar)   <- xname
      xsd           <- as.vector(xsd)
      names(xsd)    <- xname
      info$scale    <- list(mean=xbar, sd=xsd)
    }
  } else lp <- rep(kof[kmid], n)

  # Add second derivative of penalty function if needed, on the original scale
  if(! inclpen && penpres)
    info$b <- info$b - penalty.matrix

  names(kof) <- name

  stats <- NULL
  if(compstats) {
    if(p == 0) {llnull <- loglik[length(loglik)]; model.lr <- 0e0}
    else {
      llnull   <- loglik[length(loglik) - 1L]
      model.lr <- llnull - loglik[length(loglik)]
    }
    model.df <- p
    if(initial.there || maxit == 1)
      model.p <- score <- score.p <- NA
    else {
      score <- z$score
      if(model.df == 0)
        model.p <- score.p <- 1.
      else {
        model.p <- 1. - pchisq(model.lr, model.df)
        score.p <- 1. - pchisq(score,    model.df)
      }
    }

    # Effective sample size ess did not count censored observations
    # They do contain some information, especially for interval-censored ones
    # with small intervals.  Adjust ess for partial contributions.

    if(anycens) {
      nuncens <- ess
      ll      <- -2e0 * log(z$lpe)
      a       <- sum(ll[ctype == 0])
      # Compute multiplier that makes ll for uncensored obs sum to the number uncensored
      b       <- nuncens / sum(ll[ctype == 0])
      # Compute this scaled contribution of censored obs
      essc    <- b * sum(ll[ctype > 0])
      ess     <- ess + essc
    }

    r2     <- 1. - exp(- model.lr / Nw)
    r2.max <- 1. - exp(- llnull / Nw)
    r2     <- r2 / r2.max
    r2m    <- R2Measures(model.lr, model.df, Nw, ess)
    if(k > 1L) attr(lp, 'intercepts') <- kmid
    g  <- GiniMd(lp)
    ## compute average |difference| between 0.5 and the condition
    ## probability of being >= marginal median
    cump <- eval(fam[1])
    pdm <- mean(abs(cump(lp) - 0.5))
    rho <- if(anycens) NA else if(p == 0 || diff(range(lp)) == 0e0) 0e0 else cor(rank(lp), rank(y))
    ## Somewhat faster:
    ## rho <- .Fortran('rcorr', cbind(lp, y), as.integer(n), 2L, 2L, r=double(4),
    ##                 integer(4), double(n), double(n), double(n), double(n),
    ##                 double(n), integer(n), PACKAGE='Hmisc')$r[2]

    stats <- c(n, ess, length(numy), mediany, z$dmax, model.lr, model.df,
              model.p, score, score.p, rho, r2, r2m, g, exp(g), pdm)

    nam <- c("Obs", "ESS", "Distinct Y", "Median Y", "Max Deriv",
            "Model L.R.", "d.f.", "P", "Score", "Score P",
            "rho", "R2", names(r2m), "g", "gr", "pdm")
    names(stats) <- nam
    }

  info$iname <- iname
  info$xname <- xname

  retlist <- list(call              = cal,
                  freq              = numy,
                  yunique           = ylevels,
                  yupper            = yupper,             # NULL if no censoring that produced yupper > ylevels
                  Ncens             = if(isOcens) Ncens,
                  # n.risk            = if(any(ctype > 0)) n.risk,
                  yrange            = orange,             # range attribute from Ocens
                  stats             = stats,
                  coefficients      = kof,
                  var               = NULL,
                  u                 = z$u,
                  lpe               = z$lpe,
                  iter              = z$iter,
                  family            = family,
                  famfunctions      = fam,
                  deviance          = loglik,
                  non.slopes        = k,
                  interceptRef      = kmid,
                  linear.predictors = lp,
                  penalty.matrix    = if(penpres) penalty.matrix,
                  weights           = if(wtpres) weights,
                  xbar              = if(p > 0 && scale) xbar,
                  xsd               = if(p > 0 && scale) xsd,
                  info.matrix       = info,
                  fail              = FALSE)

  class(retlist) <- 'orm'
  retlist
}

ormfit <-
  function(x, y, y2, k, intcens, link, initial,
           offset=rep(0., n), wt=rep(1., n), penmat=matrix(0., p, p), opt_method='NR',
           maxit=30L, objtol=5e-4, gradtol=1e-3, paramtol=1e10, tolsolve=.Machine$double.eps,
           minstepsize=1e-2, trace=FALSE, iname, xname) {

# y  =  -1 for left censored observation
# y2 = k+1 for right censored observation
# Uncensored observations have y = y2 = 0, 1, ..., k
# intcens = 1 if there are any interval censored observations

n <- length(y)
p <- length(initial) - k

if(k > 1 && any(diff(initial[1:k]) >= 0))
  stop('initial values for intercepts are not in descending order')

storage.mode(x)       <- 'double'
storage.mode(y)       <- 'integer'
storage.mode(y2)      <- 'integer'
storage.mode(k)       <- 'integer'
storage.mode(intcens) <- 'integer'
storage.mode(p)       <- 'integer'
storage.mode(initial) <- 'double'
storage.mode(offset)  <- 'double'
storage.mode(wt)      <- 'double'
storage.mode(penmat)  <- 'double'
storage.mode(link)    <- 'integer'

rfort <- function(theta, what=3L, debug=as.integer(getOption('orm.fit.debug', 0L))) {
  p <- as.integer(length(theta) - k)
  nai <- as.integer(if(intcens) 1000000 else 0)
  if(debug) {
    a <- llist(n, k, p, y, y2, link, intcens, nai, what, debug)
    s <- sapply(a, storage.mode)
    if(any(s != 'integer')) stop(s)
    ac <- llist(x, offset, wt, penmat, theta[1:k], theta[-(1:k)],
                logL=numeric(1))
    sc <- sapply(ac, storage.mode)
    if(any(sc != 'double')) stop(s)
    g <- function(x) if(is.matrix(x)) paste(dim(x), collapse='x') else length(x)
    print(sapply(c(a, ac), g), quote=FALSE)
  }

  w <- .Fortran(F_ormll, n, k, p, x, y, y2, offset, wt, penmat,
                link=link, theta[1:k], theta[-(1:k)],
                logL=numeric(1), grad=numeric(k + p), lpe=numeric(n),
                a=matrix(0e0, (1 - intcens) * k, 2), b=matrix(0e0, p, p), ab=matrix(0e0, k, p),
                intcens, row=integer(nai), col=integer(nai), ai=numeric(nai), nai=nai, ne=integer(1),
                what=what, debug=as.integer(debug), 1L, salloc=integer(1))
  if(debug) prn(w$salloc)

  if(w$salloc == 999) {   # zero or negative probability in likelihood calculation
    w$logL = Inf          # triggers step-halving in MLE updating loop
    return(w)
  }
  if(w$salloc == 998)
    stop('Censoring values encountered that are not handled')
  if(w$salloc == 997)
    stop('More than 1,000,000 elements needed in the intercepts part of the information matrix due\n',
         'to the variety of interval censored observations')
  if(w$salloc != 0)
    stop('Failed dynamic array allocation in Fortran subroutine ormll: code ', w$salloc)
  if(intcens && what == 3L) {
    ne    <- w$ne
    w$a   <- list(row = w$row[1 : ne], col = w$col[1 : ne], a = w$ai[1 : ne])
    w$row <- w$col <- w$ai <- w$nai <- w$ne <- NULL
    }
  w
  }

  if(missing(x) || ! length(x) || p == 0) {
    x <- 0.
    p <- 0L
  }

  nv <- k + p
  if(length(initial) < nv)
    initial <- c(initial, rep(0., nv - length(initial)))
  if(trace > 2) prn(initial)

m <- function(x) max(abs(x))

if(maxit == 1) {
  w      <- rfort(initial)
  # Information matrix is negative Hessian on LL scale
  if(intcens) w$a$a <- - w$a$a else w$a <- - w$a
  info <- list(a = w$a, b = - w$b, ab = - w$ab,
               iname=iname, xname=xname)

  res <- list(coefficients = initial,
              loglik       = w$logL,
              info         = info,
              u            = w$grad,
              lpe          = w$lpe,
              dmax=m(w$grad), score=NA,
              iter=1, fail=FALSE, class='orm')
    return(res)
  }

  theta      <- initial # Initialize the parameter vector
  oldobj     <- 1e12
  score.test <- NA

  gradtol <- gradtol * n / 1e3

  # Newton-Raphson MLE with step-halving, initial draft generated by ChatGPT
  if(opt_method == 'NR') {
    for (iter in 1:maxit) {
      w <- rfort(theta)
      if(iter == 1) objf <- w$logL
      gradient <- w$grad
      hess     <- infoMxop(w[c('a', 'b', 'ab')])

      # Newton-Raphson step

      delta <- try(Matrix::solve(hess, gradient, tol=tolsolve))
      # Runs amazingly slow if Matrix:: is omitted; prob. not using Matrix
      if(inherits(delta, 'try-error')) {
        message('singular Hessian matrix')
        return(list(fail=TRUE))
      }

      if(trace > 0)
        cat('Iteration:', iter, '  -2LL:', format(objf, nsmall=4),
            '  Max |gradient|:', m(gradient),
            '  Max |change in parameters|:', m(delta), '\n', sep='')

      if(opt_method == 'NR' && is.na(score.test) && p > 0 &&
        all(theta[- (1 : k)] == 0.))
        score.test <- - gradient %*% delta

      step_size <- 1.0           # Initialize step size for step-halving

      # Step-halving loop
      while (TRUE) {
        new_theta <- theta - step_size * delta # Update parameter vector
        wd <- which(diff(new_theta[1 : k]) >= 0e0)
        if(length(wd)) {
          if(trace > 0) cat('new_theta out of order for intercepts',
              paste(wd + 1, collapse=' '), 'forced step-halving\n')
          objfnew <- Inf
        }
        else objfnew   <- rfort(new_theta, what=1L)$logL
        if(trace > 1)
          cat('Old, new, old - new -2 LL:', objf, objfnew, objf - objfnew, '\n')
        if (! is.finite(objfnew) || objfnew > objf + objtol / 100.) {
          # Objective function failed to be reduced or is infinite
          step_size <- step_size / 2e0         # Reduce the step size
          if(trace > 0) cat('Step size reduced to', step_size, '\n')
          if(step_size < minstepsize) {
            message('Step size ', step_size, ' has reduced below minstepsize=',
                    minstepsize,
                    ' without improving log likelihood; fitting stopped')
            return(list(fail=TRUE))
          }
        } else {
          theta  <- new_theta                   # Accept the new parameter vector
          oldobj <- objf
          objf   <- objfnew
          if(trace > 2) prn(theta)
          break
        }
      }

      # Convergence check - must meet 3 criteria
      if((objf <= oldobj + objtol / 100. && (oldobj - objf < objtol)) &&
        (m(gradient) < gradtol) &&
        (m(delta)    < paramtol)) {
        # Compute final information matrix (in 3 parts) since not computed
        # since Newton-Raphson updating
        w <- rfort(theta)
        if(intcens) w$a$a <- - w$a$a else w$a <- - w$a
        info <- list(a = w$a, b = - w$b, ab = - w$ab,
                    iname=iname, xname=xname)

        return(list(coef           = theta,
                    loglik         = w$logL,
                    u              = w$grad,
                    lpe            = w$lpe,
                    info           = info,
                    objchange      = oldobj - w$logL,
                    dmax           = m(w$grad),
                    maxparamchange = m(delta),
                    score          = score.test,
                    iter           = iter,
                    fail           = FALSE) )
        }
    }

    msg <- paste('Reached', maxit, 'iterations without convergence\nChange in -2LL:',
      oldobj -objf, ' Max |gradient|:', m(gradient),
      ' Max |change in parameters|:', m(delta))
    message(msg)
    return(list(fail=TRUE))

  } else {    # L-M

  lambda   <- 1e-3    # hard-wired for L-M
  oldobj   <- 1e12
  objf     <- NA      # needed in case no H_damped is ever positive definite
  w        <- rfort(theta)
  gradient <- w$grad
  H        <- infoMxop(w[c('a', 'b', 'ab')])

  for (iter in 1:maxit) {
    H_damped <- H + lambda * Matrix::Diagonal(x = Matrix::diag(H))
    delta    <- try(Matrix::solve(H_damped, gradient, tol=tolsolve))
    if(inherits(delta, 'try-error')) {
      # Increase lambda if Hessian is ill-conditioned
      lambda <- lambda * 10.
      next
      }

    theta_new <- theta - delta
    objf      <- rfort(theta_new, what=1L)$logL
    if(trace > 0)
      cat('Iteration:', iter, '  -2LL:', format(objf, nsmall=4),
          '  Max |gradient|:', m(gradient),
          '  Max |change in parameters|:', m(delta), '\n', sep='')
    if(trace > 1)
      cat('Old, new, old - new -2 LL:', oldobj, objf, oldobj - objf, '\n')

    if(is.finite(objf) &&
       (objf <= oldobj + objtol / 100. && (oldobj - objf < objtol)) &&
       (m(gradient) < gradtol) &&
       (m(delta)    < paramtol)) break

    if(is.finite(objf) && (objf < oldobj)) {
      # Accept the step and decrease lambda
      theta    <- theta_new
      oldobj   <- objf
      w        <- rfort(theta)
      gradient <- w$grad
      H        <- infoMxop(w[c('a', 'b', 'ab')])
      lambda   <- lambda / 10.
    } else {
      # Reject the step and increase lambda
      lambda <- lambda * 10.
    }
  }
  if(iter == maxit) {
    msg <- paste('Reached', maxit, 'iterations without convergence\n-2LL:',
      objf, ' Max |gradient|:', m(gradient))
    message(msg)
    return(list(fail=TRUE))
  }
  w     <- rfort(theta)
  if(intcens) w$a$a <- - w$a$a else w$a <- - w$a
  info  <- list(a = w$a, b = - w$b, ab = - w$ab,
                iname=iname, xname=xname)

  return(list(coef           = theta,
              loglik         = w$logL,
              u              = w$grad,
              lpe            = w$lpe,
              info           = info,
              objchange      = objf - w$logL,
              dmax           = m(w$grad),
              maxparamchange = m(delta),
              score          = NA,
              iter           = iter,
              fail           = FALSE) )
  }   # End M-L
}

## Note: deriv and deriv2 below are no longer used as are hard-coded into ormll
## Expressions are used because if using a function that calls plogis(),
## the .C call for plogis will can result in R losing the environment
## of the C code.  
## The 5 expression elements are cumprob, inverse, deriv, deriv2, and 
## deriv as a function only of x

## Extreme value type I dist = Gumbel maximum = exp(-exp(-x)) = MASS:::pgumbel
## Gumbel minimum = 1 - exp(-exp(x))
probabilityFamilies <-
  list(logistic =
         expression(function(x) plogis(x),
                    function(x) qlogis(x),
                    function(x, f) f*(1-f),
                    function(x, f, deriv) f*(1-3*f+2*f*f),
                    function(x) {f <- plogis(x); f*(1-f)}),
       probit =
         expression(function(x) pnorm(x),
                    function(x) qnorm(x),
                    function(x, f) dnorm(x),
                    function(x, f, deriv) - deriv * x,
                    function(x) dnorm(x)),
       loglog =
         expression(function(x) exp(-exp(-x)),
                    function(x) -log(-log(x)),
                    function(x, f) exp(-x-exp(-x)),
                    function(x, f, deriv)
                      ifelse(abs(x) > 200, 0,
                             exp(-x - exp(-x)) * (-1 + exp(-x))),
                    function(x) exp(-x-exp(-x))),
       cloglog =
         expression(function(x) 1-exp(-exp(x)),
                    function(x) log(-log(1-x)),
                    function(x, f) exp(x-exp(x)),
                    function(x, f, deriv)
                      ifelse(abs(x) > 200, 0, deriv * ( 1 - exp( x))) ,
                    function(x) exp(x-exp(x))),
       cauchit =
         expression(function(x) pcauchy(x),
                    function(x) qcauchy(x),
                    function(x, f) dcauchy(x), 
                    function(x, f, deriv) -2 * x * ((1 + x*x)^(-2)) / pi,
                    function(x) dcauchy(x))
  )

## Check:
## P(x) = plogis(x); P'(x) = P(x) - P(x)^2
## d <- function(x) plogis(x) - 3*plogis(x)^2 + 2*plogis(x)^3
## x <- seq(-3, 3, length=150)
## plot(x, d(x), type='l')
## ad <- c(NA,diff(dlogis(x))/(x[2]-x[1]))
## lines(x, ad, col='red')

