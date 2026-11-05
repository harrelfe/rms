#' Validate Predicted Probabilities Against Observed Survival Times
#'
#' @description
#' The \code{val.surv} function is useful for validating predicted survival
#' probabilities against right-censored failure times.  If \code{u} is
#' specified, one of three smoothers (see \code{method}) is used to relate
#' predicted survival probability at time \code{u} to observed survival
#' times (and censoring indicators), estimating the actual survival
#' probability at time \code{u} as a function of the estimated survival
#' probability at that time, \code{est.surv}.  If \code{est.surv} is not
#' given, \code{fit} must be specified and the \code{survest} function is
#' used to obtain the predicted values (using \code{newdata} if it is given,
#' or using the stored linear predictor values if not).  The smoother
#' (\code{hare}, \code{movStats} when \code{method="smoothkm"}, or
#' \code{loess.surv} when \code{method="loess"}) is given the sole predictor
#' \code{fun(est.surv)} where \code{fun} is given by the user or is inferred
#' from \code{fit}.  \code{fun} is the function of predicted survival
#' probabilities that one expects to create a linear relationship with the
#' linear predictors.
#'
#' @details
#' \code{method="loess"} (recommended) uses \code{\link{loess.surv}}, a
#' locally weighted smoother for right-censored data (van Klaveren et al.,
#' 2026).  Event indicators are inflated by inverse probability of
#' censoring weights obtained from a Cox model for the censoring
#' distribution, and the weighted indicators are smoothed against
#' \code{fun(est.surv)}.  By default the \code{lowess} function with
#' iterations turned off is used as the smoother, which is fast; arguments
#' such as \code{span}, \code{formula.cens}, and \code{smoother} may be
#' passed through \code{...}.
#'
#' \code{method="hare"} uses an adaptive procedure to find a linear spline of
#' \code{fun(est.surv)} in a model where the log hazard is a linear spline
#' in time \eqn{t}, and cross-products between the two splines are allowed so
#' as to not assume proportional hazards.  Thus \code{hare} assumes that the
#' covariate and time functions are smooth but not much else, if the number
#' of events in the dataset is large enough for obtaining a reliable
#' flexible fit.  Or specify \code{method="smoothkm"} to use the
#' \code{Hmisc} \code{movStats} function to compute smoothed (by default
#' using \code{supsmu}) moving window Kaplan-Meier estimates.  This method
#' is more flexible than \code{hare}.
#'
#' There are special \code{print} and \code{plot} methods when \code{u} is
#' given.  In this case, \code{val.surv} returns an object of class
#' \code{"val.survh"}, otherwise it returns an object of class
#' \code{"val.surv"}.
#'
#' If \code{u} is not specified, \code{val.surv} uses Cox-Snell (1968)
#' residuals on the cumulative probability scale to check on the calibration
#' of a survival model against right-censored failure time data.  If the
#' predicted survival probability at time \eqn{t} for a subject having
#' predictors \eqn{X} is \eqn{S(t|X)}, this method is based on the fact that
#' the predicted probability of failure before time \eqn{t}, \eqn{1 - S(t|X)},
#' when evaluated at the subject's actual survival time \eqn{T}, has a
#' uniform (0,1) distribution.  The quantity \eqn{1 - S(T|X)} is
#' right-censored when \eqn{T} is.  By getting one minus the Kaplan-Meier
#' estimate of the distribution of \eqn{1 - S(T|X)} and plotting against the
#' 45 degree line we can check for calibration accuracy.  A more stringent
#' assessment can be obtained by stratifying this analysis by an important
#' predictor variable.  The theoretical uniform distribution is only an
#' approximation when the survival probabilities are estimates and not
#' population values.
#'
#' When \code{censor} is specified to \code{val.surv}, a different
#' validation is done that is more stringent but that only uses the
#' uncensored failure times.  This method is used for type I censoring when
#' the theoretical censoring times are known for subjects having uncensored
#' failure times.  Let \eqn{T}, \eqn{C}, and \eqn{F} denote respectively
#' the failure time, censoring time, and cumulative failure time
#' distribution (\eqn{1 - S}).  The expected value of \eqn{F(T | X)} is 0.5
#' when \eqn{T} represents the subject's actual failure time.  The expected
#' value for an uncensored time is the expected value of
#' \eqn{F(T | T \leq C, X) = 0.5 F(C | X)}.  A smooth plot of
#' \eqn{F(T|X) - 0.5 F(C|X)} for uncensored \eqn{T} should be a flat line
#' through \eqn{y=0} if the model is well calibrated.  A smooth plot of
#' \eqn{2F(T|X)/F(C|X)} for uncensored \eqn{T} should be a flat line through
#' \eqn{y=1.0}. The smooth plot is obtained by smoothing the (linear
#' predictor, difference or ratio) pairs.
#'
#' Note that the Cox-Snell residual plot is not very sensitive to model lack
#' of fit.
#'
#' @param fit a fit object created by \code{cph} or \code{psm}
#' @param newdata a data frame for which \code{val.surv} should obtain
#'   predicted survival probabilities.  If omitted, survival estimates are
#'   made for all of the subjects used in \code{fit}.
#' @param S an \code{\link[survival]{Surv}} object or an \code{\link{Ocens}}
#'   object
#' @param est.surv a vector of estimated survival probabilities corresponding
#'   to times in the first column of \code{S}.
#' @param method applies if \code{u} is specified.  Use \code{"loess"} to
#'   estimate the calibration curve using \code{\link[Hmisc]{loess.surv}}
#'   (recommended; a warning is issued for other methods), \code{"hare"}
#'   (the default) to use \code{\link[polspline]{hare}}, or \code{"smoothkm"}
#'   to use \code{movStats}.
#' @param censor a vector of censoring times.  Only the censoring times for
#'   uncensored observations are used.
#' @param u a single numeric follow-up time
#' @param fun a function that transforms survival probabilities into the
#'   scale of the linear predictor.  If \code{fit} is given, and represents
#'   either a Cox, Weibull, or exponential fit, \code{fun} is automatically
#'   set to log(-log(p)).  When \code{fit} is not given and \code{method='loess'} \code{fun} is set to the identity function.
#' @param lim a 2-vector specifying limits of predicted survival
#'   probabilities for obtaining estimated actual probabilities at time
#'   \code{u}.  Default for \code{val.surv} is the limits for predictions from
#'   \code{datadist}, which for large \eqn{n} is the 10th smallest and 10th
#'   largest predicted survival probability.  For \code{plot.val.survh}, the
#'   default for \code{lim} is the range of the combination of predicted
#'   probabilities and calibrated actual probabilities.  \code{lim} is used
#'   for both axes of the calibration plot.
#' @param evaluate the number of evenly spaced points over the range of
#'   predicted probabilities.  This defines the points at which calibrated
#'   predictions are obtained for plotting.
#' @param pred a vector of points at which to evaluate predicted
#'   probabilities, overriding \code{lim}
#' @param maxdim see \code{\link[polspline]{hare}}; used only when
#'   \code{method="hare"}
#' @param x result of \code{val.surv}
#' @param xlab x-axis label.  For \code{plot.val.survh}, defaults for
#'   \code{xlab} and \code{ylab} come from \code{u} and the units of
#'   measurement for the raw survival times.
#' @param ylab y-axis label
#' @param riskdist set to \code{FALSE} to not call \code{scat1d} to draw the
#'   distribution of predicted (uncalibrated) probabilities
#' @param add set to \code{TRUE} if adding to an existing plot
#' @param scat1d.opts a \code{list} of options to pass to \code{scat1d}.
#'   By default, the option \code{nhistSpike=200} is passed so that a spike
#'   histogram is used if the sample size exceeds 200.
#' @param ... When \code{u} is given to \code{val.surv}, \code{...} represents
#'   optional arguments to \code{loess.surv} (e.g., \code{span},
#'   \code{formula.cens}, \code{smoother}), \code{hare}, or \code{movStats}.
#'   It can represent arguments to pass to \code{plot} or \code{lines} for
#'   \code{plot.val.survh}.  Otherwise, \code{...} contains optional
#'   arguments for \code{plsmo} or \code{plot}.  For \code{print.val.survh},
#'   \code{...} is ignored.
#' @param group a grouping variable.  If numeric this variable is grouped into
#'   \code{g.group} quantile groups (default is quartiles).  \code{group},
#'   \code{g.group}, \code{what}, and \code{type} apply when \code{u} is not
#'   given.
#' @param g.group number of quantile groups to use when \code{group} is given
#'   and variable is numeric.
#' @param what the quantity to plot when \code{censor} was in effect.  The
#'   default is to show the difference between cumulative probabilities and
#'   their expectation given the censoring time.  Set \code{what="ratio"} to
#'   show the ratio instead.
#' @param type Set to the default (\code{"l"}) to plot the trend line only,
#'   \code{"b"} to plot both individual subjects ratios and trend lines, or
#'   \code{"p"} to plot only points.
#' @param xlim,ylim axis limits for \code{plot.val.surv} when the
#'   \code{censor} variable was used.
#' @param datadensity By default, \code{plot.val.surv} will show the data
#'   density on each curve that is created as a result of \code{censor} being
#'   present.  Set \code{datadensity=FALSE} to suppress these tick marks drawn
#'   by \code{scat1d}.
#'
#' @return a list of class \code{"val.surv"} or \code{"val.survh"}.  Some
#'   \code{plot} methods return a \code{ggplot2} object.
#'
#' @author Frank Harrell\cr
#' Department of Biostatistics, Vanderbilt University\cr
#' fh@fharrell.com
#'
#' @references
#' van Klaveren D, Austin PC, Serruys PW, van Calster B, Harrell FE (2026).
#' Locally estimated scatterplot smoothing for censored survival data.
#' \emph{Statistics in Medicine}. \doi{10.1002/sim.70783}
#'
#' Cox DR, Snell EJ (1968): A general definition of residuals (with
#' discussion).  JRSSB 30:248--275.
#'
#' Kooperberg C, Stone C, Truong Y (1995): Hazard regression.  JASA 90:78--94.
#'
#' May M, Royston P, Egger M, Justice AC, Sterne JAC (2004): Development and
#' validation of a prognostic model for survival time data: application to
#' prognosis of HIV positive patients treated with antiretroviral therapy.
#' Stat in Med 23:2375--2398.
#'
#' Stallard N (2009): Simple tests for the external validation of mortality
#' prediction scores.  Stat in Med 28:377--388.
#'
#' @seealso \code{\link[Hmisc]{loess.surv}}, \code{\link{validate}},
#' \code{\link{calibrate}}, \code{\link[polspline]{hare}},
#' \code{\link[Hmisc]{scat1d}}, \code{\link{cph}}, \code{\link{psm}},
#' \code{\link{groupkm}}
#'
#' @examples
#' # Generate failure times from an exponential distribution
#' require(survival)
#' set.seed(123)              # so can reproduce results
#' n <- 1000
#' age <- 50 + 12*rnorm(n)
#' sex <- factor(sample(c('Male','Female'), n, rep=TRUE, prob=c(.6, .4)))
#' cens <- 15*runif(n)
#' h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
#' t <- -log(runif(n))/h
#' units(t) <- 'Year'
#' label(t) <- 'Time to Event'
#' ev <- ifelse(t <= cens, 1, 0)
#' t <- pmin(t, cens)
#' S <- Surv(t, ev)
#'
#' # First validate true model used to generate data
#'
#' # Make a smooth calibration plot for 1-year survival probability where
#' # we predict 1-year survival using the known true population survival
#' # probability.  method='loess' is the recommended method.
#' # In addition, use groupkm to show that grouping predictions into
#' # intervals and computing Kaplan-Meier estimates is not as accurate.
#'
#' s1 <- exp(-h*1)
#' w <- val.surv(est.surv=s1, S=S, u=1, method='loess')
#' plot(w, lim=c(.85,1), scat1d.opts=list(nhistSpike=200, side=1))
#' groupkm(s1, S, m=100, u=1, pl=TRUE, add=TRUE)
#'
#' # Same using hare (if polspline is available); this issues a warning
#' w <- val.surv(est.surv=s1, S=S, u=1,
#'               fun=function(p)log(-log(p)), method='hare')
#' plot(w, lim=c(.85,1))
#'
#' # Now validate the true model using residuals
#'
#' w <- val.surv(est.surv=exp(-h*t), S=S)
#' plot(w)
#' plot(w, group=sex)  # stratify by sex
#'
#'
#' # Now fit an exponential model and validate
#' # Note this is not really a validation as we're using the
#' # training data here
#' f <- psm(S ~ age + sex, dist='exponential', y=TRUE)
#' w <- val.surv(f)
#' plot(w, group=sex)
#'
#'
#' # We know the censoring time on every subject, so we can
#' # compare the predicted Pr[T <= observed T | T>c, X] to
#' # its expectation 0.5 Pr[T <= C | X] where C = censoring time
#' # We plot a ratio that should equal one
#' w <- val.surv(f, censor=cens)
#' plot(w)
#' plot(w, group=age, g=3)   # stratify by tertile of age
#'
#' @keywords models regression smooth survival
#' @concept model validation
#' @concept predictive accuracy
#' @export
val.surv <- function(fit, newdata, S, est.surv, 
                     method=c('hare', 'smoothkm', 'loess'),
                     censor, u, fun, lim, evaluate=100, pred, maxdim=5, 
                     ...)
{
  method  <- match.arg(method)
  if(method != 'loess') warning("method='loess' is recommended")
  isorm   <- ! missing(fit) && inherits(fit, 'orm')
  if(! missing(u)) {
    if(method == 'loess') fun <- function(p) p
    else if(missing(fun)) {
      if(missing(fit))
        stop('must specify fit if u is specified and fun is not')
      fun <- if(inherits(fit, 'coxph')) function(p) log(-log(p))
        else if(isorm) eval(fit$famfunctions[2])
        else if(length(it <- survreg.distributions[[fit$dist]]$quantile)) it
        else if(length(it <- survreg.auxinfo[[fit$dist]]$quantile)) it
        else stop('fun is not given and it cannot be inferred')
    }
  }
  if(missing(S)) {
    S           <- fit[['y']]
    units       <- fit$units
    if(! length(S)) stop('when S is omitted you must use y=TRUE in the fit')
  }
  else units <- attr(S, 'units')

  if(NCOL(S) == 1) S <- Surv(S)
  if(! inherits(S, c('Surv', 'Ocens'))) stop('S must be a Surv or Ocens object')
  if(inherits(S, 'Ocens')) S <- Ocens2Surv(S)
  if(ncol(S) != 2) stop('S must be a right-censored Surv object')
  
  if(missing(est.surv))
    est.surv <- if(! missing(u)) {
      if(missing(newdata))
        survest(fit, times=u, conf.int=0)$surv
      else
        survest(fit, newdata, times=u, conf.int=0)$surv
    }
    else {
      if(missing(newdata))
        survest(fit, times=S[,1], what='parallel', conf.int=0)
      else
        survest(fit, newdata, times=S[,1], what='parallel', conf.int=0)
    }
  est.surv <- unname(est.surv)
  if(! missing(u)) {
    i <- ! is.na(est.surv + S[,1] + S[,2])
    est.surv   <- est.surv[i]
    S          <- S[i,, drop=FALSE]
    curtail    <- function(x) pmin(.9999, pmax(x, .0001))
    func       <- function(x) fun(curtail(x))
    xx         <- func(est.surv)

    f <- switch(method,
                hare     = polspline::hare(S[,1], S[,2], xx,
                                           maxdim=maxdim, ...),
                smoothkm = movStats(S ~ xx, times=u, melt=TRUE, ...),
                loess    = NULL)
    if(missing(pred)) {
      if(missing(lim))
        lim <-
          datadist(est.surv)$limits[c('Low:prediction', 'High:prediction'),]
      pseq <- seq(lim[1], lim[2], length.out=evaluate)
    }
    else pseq <- pred
    
    if(method == 'hare') {
      actual    <- 1 - polspline::phare(u, func(est.surv), f)
      actualseq <- 1 - polspline::phare(u, func(pseq),     f)
    } else if(method == 'smoothkm') {
      actual    <- approx(f$xx, 1 - f$incidence, xout=func(est.surv))$y
      actualseq <- approx(f$xx, 1 - f$incidence, xout=func(pseq))$y
    } else {
      ## loess.surv returns estimated survival probabilities directly.
      ## Evaluate at the observations and at pseq in a single fit.
      k         <- length(est.surv)
      f         <- loess.surv(xx, S, u,
                              xout=c(func(est.surv), func(pseq)), ...)
      actual    <- f$y[seq_len(k)]
      actualseq <- f$y[-seq_len(k)]
    }
    w <- structure(list(p=est.surv, actual=actual,
                        pseq=pseq, actualseq=actualseq,
                        u=u, fun=fun, n=nrow(S), d=sum(S[,2]),
                        units=units), class='val.survh')
    if(method == 'hare')          w$harefit  <- f
    else if(method == 'smoothkm') w$movstats <- f
    else                          w$loessfit <- f
    return(w)
  }
  
  n <- nrow(S)
  nac <- if(! missing(fit)) fit$na.action
  if(! missing(censor) && length(censor) > 1 && ! missing(fit)) {
    if(length(censor) > n && length(nac)) {
      ## Missing observations were deleted during fit
      j <- ! is.na(naresid(nac, censor))
      censor <- censor[j]
    }
    if(length(censor) != n)
      stop("length of censor does not match # rows used in fit")
  }
  
  est.surv.censor <- lp <- NULL
  if(! missing(censor)) {
    if(missing(fit))
      stop('fit must be specified when censor is specified')
    est.surv.censor <- if(missing(newdata))
      survest(fit, times=censor, what='parallel', conf.int=0)
    else
      survest(fit, newdata, times=censor, what='parallel', conf.int=0)
    if(mc <- sum(is.na(est.surv.censor)))
      warning(paste(mc, 'observations had missing survival estimates at censoring time'))
    lp <- if(missing(newdata)) predict(fit, type='lp')
    else
      predict(fit, newdata, type='lp')
  }
  
  if(length(est.surv) != n)
    stop('length of est.surv must equal number of rows in S')
  structure(list(S=S, est.surv=est.surv,
                 censor.est.surv=if(length(est.surv.censor))
                 est.surv.censor,
                 lp=if(length(lp))lp,
                 na.action=nac), class='val.surv')
}

#' @rdname val.surv
#' @export
print.val.survh <- function(x, ...) {
  cat('\nValidation of Predicted Survival at Time=', format(x$u),
      '\tn=', x$n, ', events=', x$d, '\n\n')
  if(length(x$harefit)) {
    cat('hare fit:\n\n')
    print(x$harefit)
  }
  cat('\nFunction used to transform predictions:\n')
  cat(paste(format(x$fun), collapse=' '))
  error <- abs(x$p - x$actual)
  er <- c(mean(error, na.rm=TRUE), quantile(error, .9, na.rm=TRUE))
  erf <- format(round(er, 4))
  cat('\n\nMean absolute error in predicted probabilities:',
      erf[1],'\n')
  cat('0.9 Quantile of absolute errors               :', erf[2], '\n')
  er
}

#' @rdname val.surv
#' @export
plot.val.survh <- function(x, lim, xlab, ylab,
                           riskdist=TRUE, add=FALSE,
                           scat1d.opts=list(nhistSpike=200), ...)
  {
    if(missing(lim)) lim <- range(c(x$pseq, x$actualseq), na.rm=TRUE)
    uni <- Punits(x$units, adds=FALSE, default='unit')
    if(x$u != 1) uni <- paste0(uni, 's')
    lab <- paste('Probability of Surviving ', format(x$u), ' ', uni,
                 sep='')
    if(add)
      lines(x$pseq, x$actualseq, ...)
    else
      plot(x$pseq, x$actualseq, type='l', xlim=lim, ylim=lim,
           xlab=if(missing(xlab)) paste('Predicted', lab) else xlab,
           ylab=if(missing(ylab)) paste('Actual',    lab) else ylab)
    abline(a=0, b=1, lty=2)
    if(riskdist) do.call('scat1d', c(list(x=x$p), scat1d.opts))
    invisible()
  }
    


#' @rdname val.surv
#' @export
plot.val.surv <- function(x, group, g.group=4,
                          what=c('difference','ratio'),
                          type=c('l','b','p'),
                          xlab, ylab, xlim, ylim,
                          datadensity=TRUE, ...)
{
  S               <- x$S
  est.surv        <- x$est.surv
  censor.est.surv <- x$censor.est.surv
  what            <- match.arg(what)
  type            <- match.arg(type)
  
  n   <- length(est.surv)
  nac <- x$na.action
  
  if(! missing(group)) {
      if(length(group) > n && length(nac)) {
        ## Missing observations were deleted during fit
        j <- ! is.na(naresid(nac, est.surv))
        group <- group[j]
      }
      if(length(group) != n)
        stop("length of group does not match # rows used in fit")
      if(! is.factor(group)) group <- 
        if(is.logical(group) || is.character(group)) 
          as.factor(group) else cut2(group, g=g.group)
    }
  
  if(length(censor.est.surv)) {
    if(missing(group)) group <- rep(1, length(censor.est.surv))
    i <- S[, 2]==1
    group <- group[i]
    if(sum(i) < 2) stop('fewer than 2 uncensored observations')
    y <- switch(what,
                difference = 1-est.surv - .5*(1-censor.est.surv),
                ratio      = (1 - est.surv) / (.5 * (1 - censor.est.surv)))
    meanF <- tapply(1 - est.surv[i], group, mean, na.rm=TRUE)
    meanE <- tapply(.5 * (1 - censor.est.surv[i]), group, mean, na.rm=TRUE)
    res <- matrix(cbind(meanF, meanE), ncol=2,
                  dimnames=list(levels(group),
                    c('Mean F(T|T<C,X)', 'Expected')))
    lp <- x$lp
    lp <- lp[i]; y <- y[i]
    if(missing(xlab)) xlab <- 'Linear Predictor'
    if(missing(ylab)) ylab <-
      switch(what,
             difference='F(T|X,T<=C) - .5F(C|X)',
             ratio='F(T|X,T<=C)/.5F(C|X), C=Censoring Time')
    if(missing(xlim)) xlim <- range(lp, na.rm=TRUE)
    if(missing(ylim)) ylim <- range(y,  na.rm=TRUE)
    if(type=='l') plsmo(lp, y, group=group, xlab=xlab, ylab=ylab,
         xlim=xlim, ylim=ylim,
         datadensity=datadensity, trim=0, ...)
    else {
      plot(lp, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
      if(type!='p') plsmo(lp, y, group=group, add=TRUE,
           xlab=xlab, ylab=ylab, datadensity=datadensity, trim=0, ...)
    }
    abline(h=if(what=='difference')0 else 1, lty=2)
    return(res)
  }
  
  if(missing(xlab)) xlab <- 'Predicted Pr[T <= observed T]'
  if(missing(ylab)) ylab <- 'Fraction <= x'
  if(missing(xlim)) xlim <- 0:1
  if(missing(ylim)) ylim <- 0:1
  
  if(missing(group)) {
    nma <- ! is.na(est.surv + S[,2])
    est.surv <- est.surv[nma]
    S <- S[nma,, drop=FALSE]
    f <- survfitKM(factor(rep(1, length(est.surv))),
                   Surv(1 - est.surv,S[,2]),
                   se.fit = FALSE, conf.type = "none")
    tt <- c(0, f$time)
    ss <- c(1, f$surv)
    g <- ggplot(mapping=aes(x=tt, y=1 - ss)) + geom_step() +
      xlim(xlim) + ylim(ylim) + xlab(xlab) + ylab(ylab) +
      geom_abline(intercept = 0, slope = 1,
                  color = "black", linewidth=1, alpha=0.25)
    return(g)
  }
  
  nma <- ! (is.na(est.surv + S[,1] + S[,2]) | is.na(group))
  S        <- S[nma,, drop=FALSE]
  est.surv <- est.surv[nma]
  group    <- group[nma]
  
  ng <- length(lg <- levels(group))
  w <- NULL
  for(i in 1 : (ng + 1)) {
    s <- if(i==(ng + 1)) rep(TRUE,length(est.surv)) else group==lg[i]
    f <- survfitKM(factor(rep(1, sum(s))),
                   Surv(1 - est.surv[s], S[s, 2]),
                   se.fit = FALSE, conf.type = "none")
    w <- rbind(w, data.frame(group=c(lg, 'Overall')[i],
                             x=c(0, f$time), y=1 - c(1, f$surv)))
  }  
  g <- ggplot(w, aes(x, y, color=group)) + geom_step() +
        xlim(xlim) + ylim(ylim) + xlab(xlab) + ylab(ylab) +
        guides(color=guide_legend(title='')) +
        geom_abline(intercept = 0, slope = 1,
                    color = "black", linewidth=1, alpha=0.25)
  return(g)
}
