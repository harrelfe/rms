#' Logistic Model Fitter
#'
#' Fits a binary or propoortional odds ordinal logistic model for a given design matrix and response vector with no missing values in either.  Ordinary or quadratic penalized maximum likelihood estimation is used.
#'
#' `lrm.fit` implements a large number of optimization algorithms with the default being the quasi-Newtom method [stats::nlminb()].  For binary logistic regression without penalization iteratively reweighted least squares method in [stats::glm.fit()] is an option.  The -2 log likeilhood, gradient, and Hessian (negative information) matrix are computed in Fortran for speed.  By default, the `x` matrix is mean-centered and QR-factored to help in optimization when there are strong collinearities.  Parameter estimates and the covariance matrix are adjusted to the original `x` scale after fitting.  More detail and comparisons of the various optimization methods may be found [here](https://fharrell.com/post/mle/).  For ordinal regression with a large number of intercepts (distinct `y` values less one) you may want to use `optim_method='BFGS', compvar=FALSE` which does away with the need to compute the Hessian.  This will be helpful if statistical tests and confidence intervals are not being computed, or when only likelihood ratio tests are done.
#'
#' When there is complete separation (Hauck-Donner condition), i.e., the MLE of a coefficient is \eqn{\pm\infty}, and `y` is binary and there is no penalty, `glm.fit` may not converge because it does not have a convergence parameter for the deviance.  Setting `trace=1` will reveal that the -2LL is approaching zero but doesn't get there, relatively speaking.  In such cases the default of `nlminb` with `abstol=0.001` works well.
#' @title lrm.fit
#' @param x design matrix with no column for an intercept.  If a vector is transformed to a one-column matrix.
#' @param y response vector, numeric, categorical, or character.  For ordinal regression, the order of categories comes from `factor` levels, and if `y` is not a factor, from the numerical or alphabetic order of `y` values.
#' @param offset optional numeric vector containing an offset on the logit scale
#' @param initial vector of initial parameter estimates, beginning with the intercepts
#' @param opt_method optimization method, with possible values
#'   * `'nlminb'` : the default, a quasi-Newton method using [stats::nlminb()] which uses gradients and the Hessian.  This is a very fast and robust algorithm.
#'   * `'NR'` : standard Newton-Raphson iteration using the gradient and Hessian, with step-helving.  All three convergence criteria of `eps, gradtol, abstol` must be satisfied.  Relax some of these if you do not want to consider some of them in judging convergence.
#'   * `'glm.fit'` : for binary `y` without penalization only
#'   * `'nlm'` : see [stats::nlm()]; not highly recommended
#'   * `'BFGS'` :
#'   * `'L-BFGS-B'` :
#'   * `'CG'` :
#'   * `'Nelder-Mead'` : see [stats::optim()] for these 4 methods
#' @param maxit maximum number of iterations allowed, which means different things for different `opt_method`.  For `NR` it is the number of updates to parameters not counting step-halving steps.  When `maxit=1`, `initial` is assumed to contain the maximum likelihood estimates already, and those are returned as `coefficients`, along with `u`, `info.matrix` (negative Hessian) and `deviance`.  If `compvar=TRUE`, the `var` matrix is also returned, and `stats` are only computed if `compstats` is explicitly set to `TRUE` by the user.
#' @param reltol used by `BFGS`, `nlminb`, `glm.fit` to specify the convergence criteria in relative terms with regard to -2 LL, i.e., convergence is assume when one minus the fold-change falls below `reltol`
#' @param abstol used by `NR` (maximum absolute change in parameter estimates from one iteration to the next before convergence can be declared), `nlminb` (by default has no effect; see `abs.tol` argument; set to e.g. 0.001 for `nlminb` when there is complete separation)
#' @param gradtol used by `NR` (maximuma absolute gradient before convergence can be declared) and `nlm` (similar but for a scaled gradient)
#' @param factr see [stats::optim()] documentation for `L-BFGS-B`
#' @param eps difference in -2 log likelihood for declaring convergence with `opt_method='NR'`.  At present, the old `lrm.fit` approach of still declaring convergence even if the -2 LL gets worse by `eps/10` while the maximum absolute gradient is below 1e-9 is not implemented.  This handles the case where the initial estimates are actually MLEs, and prevents endless step-halving.
#' @param minstepsize used with `opt_method='NR'` to specify when to abandon step-halving
#' @param trace set to a positive integer to trace the iterative process.  Some optimization methods distinguish `trace=1` from `trace` higher than 1.
#' @param tol singularity criterion for inverting the -Hession matrix to get the covariance matrix, or for `opt_method='NR'` updates.
#' @param penalty.matrix a self-contained ready-to-use penalty matrix - see [lrm()].  It is \eqn{p \times p} where \eqn{p} is the number of columns of `x`.
#' @param weights a vector (same length as `y`) of possibly fractional case weights
#' @param normwt set to `TRUE` to scale `weights` so they sum to \eqn{n}, the length of `y`; useful for sample surveys as opposed to the default of frequency weighting
#' @param transx set to `FALSE` to prevent `x` from being centered and QR-factored to orthogonalize.  See [this](https://hbiostat.org/rmsc/mle#qr) for details.  When `maxit=1`, `transx` defaults to `FALSE`.
#' @param compvar set to `FALSE` to prevent the calculation of the variance-covariance matrix.
#' @param compstats set to `FALSE` to prevent the calculation of the vector of model statistics (defaults to `FALSE` when `maxit=1`)
#' @param inclpen set to `FALSE` to not include the penalty matrix in the Hessian when the Hessian is being computed on transformed `x`, vs. adding the penalty after back-transforming.  This should not matter.
#' @param initglm set to `TRUE` to compute starting values for an ordinal model by using `glm.fit` to fit a binary logistic model for predicting the probability that `y` exceeds or equals the median of `y`.  After fitting the binary model, the usual starting estimates for intercepts (log odds of cumulative raw proportions) are all adjusted so that the intercept corresponding to the median is the one from `glm.fit`.
#'
#' @return a list with the following elements:
#'   * `call`:  the R call to `lrm.fit`
#'   * `freq`:  vector of `y` frequencies
#'   * `sumty`:  vector of weighted `y` frequencies
#'   * `stats`:  vector with a large number of indexes and model parameters (`NULL` if `compstats=FALSE`):
#'     * `Obs`: number of observations
#'     * `Max Deriv`: maximum absolute gradiant
#'     * `Model L.R.`: overall model LR chi-square statistic
#'     * `d.f.`: degrees of freedom (number of non-intercepts)
#'     * `P`: p-value for the overall `Model L.R.` and `d.f.`
#'     * `C`: concordance probability between predicted probability and `y`
#'     * `Dxy`: Somer's Dxy rank correlation between predicted probability and `y`, = 2(C - 0.5)
#'     * `Gamma`:
#'     * `Tau-a`:
#'     * `R2`: documented [here](https://hbiostat/bib/r2.html/)
#'     * `Brier`: Brier score.  For ordinal models this is computed with respect the the median intercept.
#'     * `g`: g-index (Gini's mean difference of linear predictors)
#'     * `gr`: g-index on the odds ratio scale
#'     * `gp`: g-index on the probability scale
#'   * `fail`:  `TRUE` if any matrix inversion or failure to converge occurred
#'   * `coefficients`:
#'   * `var`:  variance-covariance matrix if `compvar=TRUE`
#'   * `u`:  gradient vector
#'   * `iter`:  number of iterations required.  For some optimization methods this is a vector.
#'   * `deviance`:  vector of deviances: intercepts-only, intercepts + offset (if `offset` is present), final model (if `x` is used)
#'   * `non.slopes`:  number of intercepts in the model
#'   * `linear.predictors`:  vector of linear predictors at the median intercept
#'   * `penalty.matrix`:  penalty matrix or `NULL`
#'   * `weights`:  `weights` or `NULL`
#'   * `xbar`:  vector of column means of `x`, or `NULL` if `transx=FALSE`
#'   * `xtrans`:  input value of `transx`
#'   * `R`: R matrix from QR to be used to rotate parameters back to original scale in the future
#'   * `Ri`: inverse of `R`
#'   * `opt_method`:  input value
#' @md
#' @author Frank Harrell <fh@fharrell.com>
#' @export
#' @keywords models regression logistic
#'
#' @examples
#' \dontrun{
#' # Fit an additive logistic model containing numeric predictors age,
#' # blood.pressure, and sex, assumed to be already properly coded and
#' # transformed
#'
#' fit <- lrm.fit(cbind(age,blood.pressure,sex), death)
#' }
#' @seealso [lrm()], [stats::glm()], [cr.setup()], [gIndex()]
lrm.fit <-
  function(x, y, offset=0, initial,
    opt_method=c('nlminb', 'NR', 'glm.fit', 'nlm', 'BFGS', 'L-BFGS-B', 'CG', 'Nelder-Mead'),
    maxit=50, reltol=1e-10,
    abstol=if(opt_method=='NR') 1e-5 else 0e0,
    gradtol=1e-5, factr=1e7, eps=5e-4,
    minstepsize=1e-4, trace=0,
    tol=1e-13, penalty.matrix=NULL, weights=NULL, normwt=FALSE,
    transx=maxit > 1, compvar=TRUE, compstats=maxit > 1,
    inclpen=TRUE, initglm=FALSE)
{

  cal                     <- match.call()
  opt_method_given        <- ! missing(opt_method)
  opt_method              <- match.arg(opt_method)
  penpres                 <- length(penalty.matrix) && any(penalty.matrix != 0)
  original.penalty.matrix <- if(penpres) penalty.matrix

  ftdb <- as.integer(getOption('lrm.fit.debug', 0L))
  if(ftdb) {
    w <- llist(opt_method, length(x), dim(x), length(y), length(offset),
               if(! missing(initial)) initial, compvar, compstats, ftdb, maxit, labels=FALSE)
    prn(w, file='/tmp/z')
  }
  db <- if(! ftdb) function(...) {} else function(...) cat(..., '\n', file='/tmp/z', append=TRUE)

  n <- length(y)

  wtpres <- TRUE
  if(! length(weights)) {
    wtpres <- FALSE
    normwt <- FALSE
    weights <- rep(1, n)
  }
  if(length(weights) != n) stop('length of wt must equal length of y')
  if(normwt) weights <- weights * n / sum(weights)
  storage.mode(weights) <- 'double'

  R <- Ri <- NULL

  initial.there <- ! missing(initial)
  if(missing(x) || length(x) == 0) {
    p     <- 0L
    xname <- NULL
    x     <- 0e0
  } else {
    if(! is.matrix(x)) x <- as.matrix(x)
    xname <- colnames(x)
    dx <- dim(x)
    p  <- dx[2]
    if(dx[1] != n) stop("x and y must have same number of rows")

    # Center X columns then QR-transform them
    # See https://betanalpha.github.io/assets/case_studies/qr_regression.html
    if(transx) {
      w    <- qrxcenter(x)
      xbar <- w$xbar
      x    <- w$x
      R    <- w$R
      Ri   <- w$Ri
      w    <- NULL
    }

    storage.mode(x) <- "double"
    if(length(xname) == 0) xname <- paste("x[", 1 : ncol(x), "]", sep="")

    # Modify penalty matrix so that it applies to the QR-transformed beta
    # With R (p x p), beta = original scale beta, gamma = beta for transformed x
    # beta = R x gamma
    # Penalty on LL ignoring -0.5 factor is beta' P beta =
    #   (R x gamma)' P (R x gamma) = gamma' R'PR gamma
    # So change P to R'PR
    if(transx && penpres)
      penalty.matrix <- t(R) %*% penalty.matrix %*% R
  }

  if(! is.factor(y)) y <- as.factor(y)
  numy    <- table(y)
  ylevels <- names(numy)
  k       <- length(ylevels) - 1
  y       <- as.integer(y) - 1  # 0, 1, 2, 3, ... k
  ymed    <- max(1, round(median(y)))

 
  if(opt_method == 'glm.fit' && (k > 1 || penpres))
    stop('opt_method="glm.fit" only applies when k=1 and there is no penalty')

  offset <- rep(offset, length=n)
  ofpres <- ! all(offset == 0)
  if(ofpres && length(offset) != n) stop("offset and y must have same length")
  storage.mode(offset) <- "double"

  if(n < 3) stop("must have >=3 non-missing observations")
  nv   <- p + k

  sumwty <- tapply(weights, y, sum)
  sumwt  <- sum(sumwty)
  if(! wtpres && any(numy != sumwty)) stop('program logic error 1')
  sumw <- if(normwt) numy else as.integer(round(sumwty))

  if(missing(initial)) {
    ncum    <- rev(cumsum(rev(sumwty)))[2 : (k + 1)]
    pp      <- ncum / sumwt
    initial <- qlogis(pp)
    if(ofpres) initial <- initial - mean(offset)
  }
  if(length(initial) < nv)
    initial <- c(initial, rep(0, nv - length(initial)))

  loglik <- -2 * sum(sumwty * logb(sumwty / sum(sumwty)))

  if(p > 0) {
    if(! penpres) penalty.matrix <- matrix(0, nrow=p, ncol=p)
    if(nrow(penalty.matrix) != p || ncol(penalty.matrix) != p)
      stop(paste("penalty.matrix does not have", p, "rows and columns"))
  }
  else
    penalty.matrix <- 0
  storage.mode(penalty.matrix) <- 'double'

  cont <- switch(opt_method,
           BFGS      = list(trace=trace, maxit=maxit, reltol=reltol),
          'L-BFGS-B' = list(trace=trace, maxit=maxit, factr=factr),
          nlminb     = list(trace=trace, iter.max=maxit, eval.max=maxit,
                            rel.tol=reltol, abs.tol=abstol, xf.tol=1e-16),
          glm.fit    = list(epsilon=reltol, maxit=maxit, trace=trace),
                       list(trace=trace, maxit=maxit) )

  info  <- NULL
  statsnull <- c(n, 0, 0, 0, NA, rep(0, 6), rep(NA, 4))

  if(p == 0 & ! ofpres) {
    loglik <- rep(loglik, 2)
    res <- list(coef=initial,
                u=rep(0, k),
                stats=statsnull)
  }

  envi <- .GlobalEnv

  ialpha = 1 : k

  logl <- function(parm, ...) {
    p     <- length(parm) - k
    ibeta <- if(p == 0) integer(0) else - (1 : k)
    db('1', n, k, p)
    .Fortran(F_lrmll, n, k, p, x, y, offset, weights, penalty.matrix,
             as.double(parm[ialpha]), as.double(parm[ibeta]),
             logL=numeric(1), numeric(nv), numeric(0), 1L, ftdb, 0L)$logL
  }

  grad <- function(parm, ...) {
    p     <- length(parm) - k
    ibeta <- if(p == 0) integer(0) else - (1 : k)
    db('2', n, k, p)
    g <- .Fortran(F_lrmll, n, k, p, x, y, offset, weights, penalty.matrix,
                  as.double(parm[ialpha]), as.double(parm[ibeta]),
                  logL=numeric(1), grad=numeric(k + p), numeric(0),
                  2L, ftdb, 0L)$grad
    # Save last computed gradient
    if(! outputs_gradient) assign('.lrm_gradient.', g, envir=envi)
    -2e0 * g
  }

  Hess <- function(parm, ...) {
    p     <- length(parm) - k
    ibeta <- if(p == 0) integer(0) else - (1 : k)
    # Last argument to lrmll = penhess = 1 so that penalty matrix is
    # included in Hessian so that it's used in Newton-Raphson updating.
    # Later for evaluating the Hessian a final time for getting the
    # covariance matrix we'll want to suppress penhess
    db('3', n, k, p)
    -2e0 * .Fortran(F_lrmll, n, k, p, x, y, offset, weights, penalty.matrix,
                    as.double(parm[ialpha]), as.double(parm[ibeta]),
                    logL=numeric(1), grad=numeric(k + p),
                    h=matrix(0e0, nrow=k + p, ncol=k + p),
                    3L, ftdb, 1L)$h
  }

  outputs_gradient <- opt_method %in% c('NR', 'nlm')

  mle <- function(p, init) {

    if(opt_method == 'NR') {
      opt <- newtonr(init, logl, grad, Hess,
                     objtol      = eps,
                     gradtol     = gradtol,
                     paramtol    = abstol,
                     minstepsize = minstepsize,
                     tolsolve    = tol,
                     maxit       = maxit,
                     trace       = trace)
      ll  <- opt$obj
      cof <- opt$param
      gr  <- -0.5 * opt$grad

      # Compute info matrix (- Hessian)
      # Hess * -0.5 = Hessian
      info <- 0.5 * Hess(cof)      
      # This is the information matrix for (delta, gamma) on the centered and QR-transformed
      # x without including the penalty matrix (which will be added below after
      # transformations) if inclpen is FALSE.

      it  <- opt$iter
      con <- opt$code
      ok  <- con == 0
    }
    else if(opt_method == 'nlm') {
      obj <- function(parm, ...) {
        ob <- logl(parm, ...)
        attr(ob, 'gradient') <- grad(parm, ...)
        attr(ob, 'hessian')  <- Hess(parm, ...)
        ob
      }
      opt  <- nlm(obj, init, gradtol=gradtol, iterlim=maxit, print.level=trace,
                 hessian=FALSE)
      ll   <- opt$minimum
      cof  <- opt$estimate
      gr   <- -0.5 * opt$gradient
      info <- 0.5 * Hess(cof)
      it   <- opt$iterations
      con  <- opt$code
      ok   <- con <= 3
    }
    else if(opt_method == 'nlminb') {
      opt  <- nlminb(start=init, objective=logl, gradient=grad, hessian=Hess,
                     control=cont)
      ll   <- opt$objective
      cof  <- opt$par
      gr   <- .lrm_gradient.
      info <- 0.5 * Hess(cof)
      it   <- c(iterations=opt$iterations, evaluations=opt$evaluations)
      con  <- opt$convergence
      ok   <- con == 0
    }
    else if(opt_method == 'glm.fit') {
      f <- try(glm.fit(cbind(1e0, x), y, family=binomial(),
                    weights=weights, offset=offset,
                    singular.ok=FALSE,
                    control=cont))
      if(inherits(f, 'try-error')) {
        message('glm.fit failed; consider using default opt_method')
        return(structure(list(fail=TRUE), class='lrm'))
      }
      ll   <- f$deviance
      cof  <- f$coefficients
      # grad() took 0.001s for n=10000 p=50
      gr   <- grad(cof)
      info <- crossprod(qr.R(f$qr))
      it   <- f$iter
      con  <- 1 - f$converged
      ok   <- con == 0
    } else {
      opt <- optim(init, logl, grad,
                   method=opt_method, control=cont, hessian=FALSE)
      ll   <- opt$value
      cof  <- opt$par
      gr   <- .lrm_gradient.
      info <- 0.5 * Hess(cof)
      it   <- opt$counts
      con  <- opt$convergence
      ok   <- con == 0
    }

    if(! ok) {
      msg  <- opt$message
      if(! length(msg))
        msg <- paste('Iterations, function and gradient evaluations:',
                      paste(it, collapse=' '))
      message('optimization did not converge for lrm.fit', '\n', msg,
              '\nconvergence code ', con)
      return(structure(list(fail=TRUE), class='lrm'))
    }
    ibeta <- if(p == 0) integer(0) else - (1 : k)
    list(logl=ll, alpha=cof[ialpha], beta=cof[ibeta], u=gr,
         info=info, iter=it)
  }

stats <- function(lp, u) {
    llnull   <- loglik[length(loglik) - 1]
    model.lr <- llnull - loglik[length(loglik)]
    model.df <- p
    model.p  <- if(model.df > 0) 1 - pchisq(model.lr, model.df) else 1e0

    r2     <- 1 - exp(- model.lr / sumwt)
    r2.max <- 1 - exp(- llnull   / sumwt)
    r2     <- r2 / r2.max
    r2m    <- R2Measures(model.lr, model.df, sumwt, sumwty)
    prob   <- plogis(lp)
    event  <- y > (ymed - 1)
    # B <- mean((prob - event)^2)
    B     <- sum(weights*(prob - event)^2) / sum(weights)
    g     <- GiniMd(lp)
    gp    <- GiniMd(prob)
    a     <- suppressWarnings(survival::concordancefit(y, lp, reverse=FALSE))
    conc  <- a$count['concordant']
    disc  <- a$count['discordant']
    tiedy <- a$count['tied.y']
    fn    <- as.double(n)
    pairs <- fn * (fn - 1e0) / 2e0
    rankstats <- c(C     = a$concordance,
                   Dxy   = (conc - disc) / (pairs - tiedy),
                   Gamma = (conc - disc) / (conc + disc),
                   Tau_a = (conc - disc) / pairs)
    stats <- c(n, max(abs(u)), model.lr, model.df,
               model.p, rankstats,
               r2, r2m, B, g, exp(g), gp)

    nam <- c("Obs", "Max Deriv",
             "Model L.R.", "d.f.", "P", "C", "Dxy",
             "Gamma", "Tau-a", "R2", names(r2m), "Brier", "g", "gr", "gp")

    names(stats) <- ifelse(nam == 'R2m', names(r2m), nam)

    if(wtpres) stats <- c(stats, 'Sum of Weights'=sumwt)
    stats
  }
  
  storage.mode(n) <- storage.mode(k) <- storage.mode(p) <- storage.mode(y) <- 'integer'
  v <- vc <- NULL

  if(maxit == 1) {
    loglik <- c(loglik, logl(initial))
    if(p == 0) lpmid <- rep(initial[ymed], n)
    else {
      lp     <- matxv(x, initial, kint=1)
      lpmid  <- lp - initial[1] + initial[ymed]
    }

    # Hess() returns Hessian on -2LL scale
    # Information matrix is negative Hessian on LL scale
    u <- -0.5 * grad(initial)
    res <- list(coefficients = initial,
                deviance     = loglik,
                info.matrix  =  0.5 * Hess(initial),
                u            = u,
                stats        = if(compstats) stats(lpmid, u),
                maxit = 1, fail=FALSE, class='lrm')
    if(compvar) {
        vc <- try(solve(res$info.matrix, tol=tol))
        if(inherits(vc, 'try-error')) {
          cat(vc)
          return(structure(list(fail=TRUE), class='lrm'))
        }
      res$var <- vc
    }
    return(res)
  }

  if(ofpres) {
    # Fit model with only intercept(s) and offset
    ml <- mle(0L, initial[1 : k])
    if(inherits(ml, 'lrm')) return(ml)
    loglik  <- c(loglik, ml$logl)
    res     <- list(coef=ml$alpha, u=ml$u, info=ml$info, iter=ml$iter)
    initial[1 : k] <- ml$alpha
    }

  # Fit full model
  if(p > 0) {
    # If k > 1, do like MASS::polr in using glm.fit to get
    # initial parameter estimates for ordinal model
    if(k > 1 && initglm) {
      f <- try(glm.fit(cbind(1e0, x), 1L*(y >= ymed), family=binomial(),
                    weights=weights, offset=offset,
                    singular.ok=FALSE))
      if(inherits(f, 'try-error') || ! f$converged) {
        message('glm.fit failed with initglm')
        return(structure(list(fail=TRUE), class='lrm'))
      }
      # Intercept in this binary Y fit corresponds to Y=ymed
      initial[1 : k] <- initial[1 : k] - initial[ymed] + f$coefficients[1]
      initial[-(1 : k)] <- f$coefficients[-1]
    }
      
      
    ml <- mle(p, initial)
    if(inherits(ml, 'lrm')) return(ml)
    loglik <- c(loglik, ml$logl)
    delta  <- ml$alpha
    gamma  <- ml$beta
    info   <- ml$info

    if(transx) {
      # Let M = matrix(1, nrow=k, ncol=1) %*% matrix(xbar, nrow=1)
      # beta  = R gamma p x p * p x 1 = p x 1
      # alpha = delta - M beta = delta - M R gamma = kx1 - kxp pxp px1
      #
      # Transform the information matrix from the (delta, gamma) scale to (alpha, beta)
      # See https://hbiostat.org/rmsc/mle#qr and https://stats.stackexchange.com/questions/657210
      # In the future see the R madness package
      M <- matrix(1, nrow=k, ncol=1) %*% matrix(xbar, nrow=1)
      J <- rbind(cbind(diag(k),                      M ),
                 cbind(matrix(0e0, nrow=p, ncol=k),  Ri))
      info <- t(J) %*% info %*% J
      }
      
    # Add second derivative of penalty function if needed, on the original scale
    if(! inclpen && penpres)
      info[-(1:k), -(1:k)] <- info[-(1:k), -(1:k)] - original.penalty.matrix

  if(compvar) {
      vc <- try(solve(info, tol=tol))
      if(inherits(vc, 'try-error')) {
        cat(vc)
        return(structure(list(fail=TRUE), class='lrm'))
      }
    }
    # Save predictions before reverting since original x not kept
    # x and kof both translated -> x * kof on original scale
    kof    <- c(delta, gamma)
    lp     <- matxv(x, kof, kint=1)
    lpmid  <- lp - kof[1] + kof[ymed]

    if(transx) {
      # Revert parameters to original x space without centering.  R is p x p
      # x is reproduced by x %*% solve(R) + matrix(1, nrow=n) %*% matrix(xbar, ncol=1)
      # If beta is p x 1 the rotation is R x beta instead of beta %*% t(R)
      delta <- matrix(delta, ncol=1)
      gamma <- matrix(gamma, ncol=1)
      beta  <- R %*% gamma
      alpha <- delta - sum(beta * xbar)
      # In matrix form alpha = delta - M beta = delta - M R gamma where
      # M = matrix(1, nrow=k, ncol=1) %*% matrix(xbar, nrow=1)
    } else {alpha <- delta; beta <- gamma}
    res <- list(coef=c(alpha, beta), u=ml$u, iter=ml$iter)
  } else {
    lp    <- rep(res$alpha[1],    n)
    lpmid <- rep(res$alpha[ymed], n)
  }

  name <- if(k == 1) 'Intercept' else paste('y>=', ylevels[2 : (k + 1)], sep="")
  name <- c(name, xname)
  kof  <- res$coef

  names(kof) <- name
  if(length(res$u)) names(res$u) <- name
  if(! is.null(info)) dimnames(info) <- list(name, name)
  if(! is.null(vc))   dimnames(vc)   <- list(name, name)
  if(length(R))  dimnames(R) <- dimnames(Ri) <- list(xname, xname)

  retlist <-
     list(call              = cal,
          freq              = numy,
          sumwty            = if(wtpres) sumwty,
          stats             = if(p == 0) statsnull else if(compstats) stats(lpmid, res$u),
          fail              = FALSE,
          coefficients      = kof,
          info.matrix       = info,
          var               = vc,
          u                 = res$u,
          iter              = res$iter,
          deviance          = loglik,
          non.slopes        = k,
          linear.predictors = lp,
          penalty.matrix    = if(p > 0 && penpres) original.penalty.matrix,
          weights           = if(wtpres) weights,
          xbar              = if(transx & p > 0) xbar,
          xtrans            = transx,
          R                 = if(transx & p > 0) R,
          Ri                = if(transx & p > 0) Ri,
          opt_method        = opt_method)

  class(retlist) <- 'lrm'
  retlist
}

# Newton-Raphson MLE with step-halving, initial draft generated by ChatGPT
# The Hessian needs to be computed by the user on the
# final param values after newtonr completes

newtonr <- function(init, obj, grad, hessian,
                    objtol = 5e-4, gradtol = 1e-5, paramtol = 1e-5,
                    minstepsize = 5e-4,
                    tolsolve=1e-7, maxit = 100, trace=0) {

  m <- function(x) max(abs(x))

  theta <- init # Initialize the parameter vector
  oldobj <- 1e10
  objf   <- obj(theta)

  for (iter in 1:maxit) {
    gradient <- grad(theta)     # Compute the gradient vector
    hess     <- hessian(theta)  # Compute the Hessian matrix

    delta <- try(qr.solve(hess, gradient, tol=tolsolve)) # Compute the Newton-Raphson step
    if(inherits(delta, 'try-error'))
      return(list(code=2, message='singular Hessian matrix'))

    if(trace > 0)
      cat('Iteration:', iter, '  -2LL:', objf, '  Max |gradient|:', m(gradient),
          '  Max |change in parameters|:', m(delta), '\n', sep='')

    step_size <- 1                # Initialize step size for step-halving

    # Step-halving loop
    while (step_size > minstepsize) {
      new_theta <- theta - step_size * delta # Update parameter vector
      objfnew <- obj(new_theta)
      if (objfnew > objf) {     # Objective function failed to be reduced
        step_size <- step_size / 2e0         # Reduce the step size
        if(trace > 0) cat('Step size reduced to', step_size, '\n')
      } else {
        theta  <- new_theta                   # Accept the new parameter vector
        oldobj <- objf
        objf   <- objfnew
        break
      }
    }

    # Convergence check - must meet 3 criteria
    if((objf <= oldobj && (oldobj - objf < objtol)) &&
       (m(gradient) < gradtol) &&
       (m(delta)    < paramtol))
        return(list(param          = theta,
                    obj            = objf,
                    grad           = gradient,
                    objchange      = oldobj - objf,
                    maxgrad        = m(gradient),
                    maxparamchange = m(delta),
                    iter=iter, code=0, message=''))
  }

  msg <- paste('Reached', maxit, 'iterations without convergence\nChange in -2LL:',
    oldobj -objf, ' Max |gradient|:', m(gradient),
    ' Max |change in parameters|:', m(delta))

  list(code = 1, message=msg)
}
