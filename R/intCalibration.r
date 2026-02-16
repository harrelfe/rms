#' Check Parallelism Assumption of Ordinal Semiparametric Models
#'
#' For all the observations used a model fit, computes the estimated probability that Y is greater than each of a number of cutoffs, and compares this to smoothed estimated probabilities as a function of predicted probabilities, to obtain internal model calibration plots with multiple cutpoints.  When Y is uncensored these are smoothed moving empirical cumulative distribution function estimates, and when Y has censored observations these are smoothing moving Kaplan-Meier estimates.  [Hmisc::movStats()] is used to do the moving overlapping window calculations.  When `hare=TRUE`, adaptive linear spline hazard regression estimates are also made, using [polspline::hare()].  When `ordsurv=TRUE`, adaptive ordinal regression estimates are made in addition.
#'
#' These plots are plots of calibration-in-the-small.  Alternate calibration-in-the-small plots may be obtained by specifying a predictor variable `x` against which to plot both predicted and observed probabilties as a function of `x`.  This is the only place in the `rms` package where the "total effect" of a predictor is estimated instead of a partial effect.  When `x` varies and moving overlapping windows of predicted and observed exceedance probabilities are estimated, if `x` is collinear with other predictors, they will "come along for the ride".
#'
#' The function also prints information on calibration-in-the-large, i.e., the mean predicted probability of being beyond each cutpoint vs. the overall proportion of observations above that cutpoint.  This is when `x` is not given.
#'
#' @param fit a fit object for which there is a [survest()] method, with `x=TRUE, y=TRUE` in effect
#' @param ycuts a vector of cutpoints on Y
#' @param m used when `ycuts` is not given.  The lowest cutoff is chosen as the first Y value having at meast `m` uncensored observations to its left, and the highest cutoff is chosen so that there are at least `m` uncensored observations to the right of it.  Cutoffs are equally spaced between these values in terms of number of uncensored observations.  If omitted, `m` is set to the minimum of 50 and one quarter of the uncensored sample size.
#' @param x a variable for which calibration-in-the-small is desired, instead of plotting predicted vs. observed probabilities.  `x` will typically be chosen by virtue of being a strong predictor (such that lack of fit will matter more) but doesn't have to be in the model.
#' @param onlydata set to `TRUE` to return a data frame suitable for plotting instead of actually plotting
#' @param eps,bass,tsmooth,hare,ordsurv see [Hmisc::movStats()]
#' @param dec number of digits to the right of the decimal place to which to round computed `ycuts`
#' @param xlab x-axis label with default constructed from the Y-variable name in the model fit (y-axis label when `x` is specified)
#' @param ylab y-axis label
#' @param nrow if `hare=TRUE` or `ordsurv=TRUE`, the number of rows in the graph (must be 1 or 2)
#' @param ... other arguments passed to [Hmisc::movStats()].  To control the number of knots for `ordsurv=TRUE` specify `k=` here.
#' @returns `ggplot2` object or a data frame
#' @export
#' @md
#' @author Frank Harrell
#' @examples
#' \dontrun{
#' getHdata(nhgh)
#' f <- orm(gh ~ rcs(age, 4), data=nhgh, family='loglog', x=TRUE, y=TRUE)
#' intCalibration(f, ycuts=c(5, 5.5, 6, 6.5))
#' f <- update(f, family='cloglog')
#' intCalibration(f, ycuts=c(5, 5.5, 6, 6.5))
#' intCalibration(f, ycuts=c(5, 6, 7), x=nhgh$age)
#' }
intCalibration <-
  function(fit, ycuts, m, x, onlydata=FALSE,
           eps=25, bass=9, tsmooth='lowess', hare=TRUE, ordsurv=TRUE,
           dec=4, xlab=bquote(hat(P)(.(yname) > y)),
           ylab='Nonparametric Estimate', nrow=1, ...) {

  Y <- fit[['y']]
  if(! length(Y)) stop('requires y=TRUE specified to fitting function')

  isocens <- inherits(Y, 'Ocens')
  if(isocens) Y <- Ocens2Surv(Y)
  else if(! survival::is.Surv(Y)) Y <- survival::Surv(Y)
  yname  <- fit$yname
  yunits <- fit$units
  if(! length(yunits)) yunits <- ''

  # Find cuts such that there are m uncensored observations beyond outer cuts and
  # between interior cuts

  if(missing(ycuts)) {
    yu <- Y[Y[, 2] == 1, 1]
    nu <- length(yu)
    if(missing(m)) m <- min(50, floor(nu / 4))
    if(nu < 2 * m) stop('number of uncensored observations ', nu,
                                ' < 2 * m =', 2 * m)
    ycuts <- cutGn(yu, m=m, what='summary')[, 'max']
    ycuts <- round(ycuts[- length(ycuts)], dec)
  }

  s <- survest(fit, times=ycuts, conf.int=0)

  if(! missing(x)) {
    xname <- deparse(substitute(x))
    vlab  <- label(x)
    if(vlab == '') xvab <- xname
    nac <- fit$na.action
    if(length(nac) && length(nac$omit)) x <- x[- nac$omit]
    if(length(x) != NROW(Y))
      stop('length of x after removing observations discarded during the fit (', length(x), ')\n',
           'is not equal to the number of observations used in the fit (', NROW(Y), ')')
    xdisc <- is.character(x) || is.factor(x) || length(unique(x)) < 10
    if(xdisc) {hare <- ordsurv <- FALSE}
    R <- NULL

    for(y in ycuts) {
      sy    <- s[s$time == y,, drop=FALSE]
      spred <- movStats(surv ~ x, data=sy, melt=TRUE, discrete=xdisc,
                        stat=function(x) list(Mean = mean(x)),
                        tunits=fit$units, tsmooth=tsmooth, hare=hare, ordsurv=ordsurv,
                        eps=eps, bass=bass, ...)
      spred$y    <- y
      spred$Type <- 'Predicted'
      sobs  <- movStats(Y ~ x, times=y, melt=TRUE, discrete=xdisc,
                        tunits=fit$units, tsmooth=tsmooth, hare=hare, ordsurv=ordsurv,
                        eps=eps, bass=bass, ...)
      sobs$y    <- y
      sobs$Type <- if(xdisc) 'Observed'
       else c(Moving='Observed (moving K-M)', HARE='Observed (HARE)', orm='Observed (orm)')[sobs$Type]
      sobs$surv <- unclass(1 - sobs$incidence)
      sobs$incidence <- NULL
      R <- rbind(R, spred, sobs)
    }
    i <- R$surv >= 0 & R$surv <= 1
    R <- R[i, ]
    if(onlydata) return(R)
    g <- ggplot(R, aes(x=.data$x, y=.data$surv, col=.data$Type)) +
           xlab(vlab) + ylab(xlab) + guides(color=guide_legend(title=''))
    if(xdisc) g <- g + geom_point() else g <- g + geom_line()
    g <- g + facet_wrap(~ paste0(.data$y, if(yunits != '') paste0('-', yunits)))
    return(g)
  }

  km_overall <- km.quick(Y, times=ycuts)
  mp <- with(s, tapply(surv, time, mean))
  d <- data.frame(y=ycuts, 'Mean Predicted P(Y > y)'=mp, 'Observed P(Y > y)'=km_overall,
                  check.names=FALSE)
  cat('\nCalibration-in-the-large:\n\n')
  print(d, digits=4, row.names=FALSE)

  R <- NULL
  for(y in ycuts) {
    sy <- s[s$time == y,,drop=FALSE]
    km <- movStats(Y ~ surv, times=y, data=sy, melt=TRUE,
                   tunits=fit$units, tsmooth=tsmooth, hare=hare, ordsurv=ordsurv,
                   eps=eps, bass=bass, ...)
    hh <- function(x) sum(! is.na(x))
    i <- km$incidence >= 0 & km$incidence <= 1
    R <- rbind(R, km[i,,drop=FALSE])
  }
  if(onlydata) return(R)
  g <- ggplot(R, aes(x=.data$surv, y=1 - .data$incidence, color=.data$Statistic)) + geom_line() +
    geom_abline(intercept=0, slope=1, alpha=0.3) +
    xlab(xlab) + ylab(ylab) + guides(color=guide_legend(title=expression(y)))
  if(hare || ordsurv) g <- g + facet_wrap(~ .data$Type, nrow=nrow)
  g
}
