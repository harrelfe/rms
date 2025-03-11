#' Check Parallelism Assumption of Ordinal Semiparametric Models
#'
#' For all the observations used a model fit, computes the estimated probability that Y is greater than each of a number of cutoffs, and compares this to smoothed estimated probabilities as a function of predicted probabilities, to obtain internal model calibration plots with multiple cutpoints.  When Y is uncensored these are smoothed moving empirical cumulative distribution function estimates, and when Y has censored observations these are smoothing moving Kaplan-Meier estimates.  [Hmisc::movStats()] is used to do the moving overlapping window calculations.  When `hare=TRUE`, adaptive linear spline hazard regression estimates are also made, using [polspline::hare()].
#' 
#' The function also prints information on calibration-in-the-large, i.e., the mean predicted probability of being beyond each cutpoint vs. the overall proportion of observations above that cutpoint.
#'
#' @param fit a fit object for which there is a [survest()] method, with `x=TRUE, y=TRUE` in effect
#' @param ycuts a vector of cutpoints on Y
#' @param m used when `ycuts` is not given.  The lowest cutoff is chosen as the first Y value having at meast `m` uncensored observations to its left, and the highest cutoff is chosen so that there are at least `m` uncensored observations to the right of it.  Cutoffs are equally spaced between these values in terms of number of uncensored observations.  If omitted, `m` is set to the minimum of 50 and one quarter of the uncensored sample size.
#' @param onlydata set to `TRUE` to return a data frame suitable for plotting instead of actually plotting
#' @param eps,bass,tsmooth,hare see [Hmisc::movStats()]
#' @param dec number of digits to the right of the decimal place to which to round computed `ycuts`
#' @param xlab x-axis label with default constructed from the Y-variable name in the model fit
#' @param ylab y-axis label
#' @param nrow if `hare=TRUE`, the number of rows in the graph (must be 1 or 2)
#' @param ... other arguments passed to [Hmisc::movStats()]
#' @returns `ggplot2` object or a data frame
#' @export
#' @md
#' @author Frank Harrell
#' @examples
#' \dontrun{
#' getHdata(nhgh)
#' f <- orm(gh ~ rcs(age, 4) + ran, data=nhgh, family='loglog', x=TRUE, y=TRUE)
#' intCalibration(f, ycuts=c(5, 5.5, 6, 6.5))
#' f <- update(f, family='cloglog')
#' intCalibration(f, ycuts=c(5, 5.5, 6, 6.5))
#' }
intCalibration <-
  function(fit, ycuts, m, onlydata=FALSE, 
           eps=25, bass=9, tsmooth='lowess', hare=TRUE,
           dec=4, xlab=bquote(hat(P)(.(yname) > y)),
           ylab='Nonparametric Estimate', nrow=1, ...) {
 
  Y <- fit[['y']]
  if(! length(Y)) stop('requires y=TRUE specified to fitting function')
  
  isocens <- inherits(Y, 'Ocens')
  if(isocens) Y <- Ocens2Surv(Y)
  else if(! survival::is.Surv(Y)) Y <- survival::Surv(Y)

  yname <- fit$yname

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

  s          <- survest(fit, times=ycuts, conf.int=0)
  km_overall <- km.quick(Y, times=ycuts)
  mp <- with(s, tapply(surv, time, mean))
  d <- data.frame(y=ycuts, 'Mean Predicted P(Y > y)'=mp, 'Observed P(Y > y)'=km_overall,
                  check.names=FALSE)
  cat('\nCalibration-in-the-large:\n\n')
  print(d, digits=4, row.names=FALSE)

  R <- NULL
  for(y in ycuts) {
    sy <- subset(s, time == y)
    km <- movStats(Y ~ surv, times=y, data=sy, melt=TRUE,
                   tunits=fit$units, tsmooth=tsmooth, hare=hare,
                   eps=eps, bass=bass, ...)
    i <- km$incidence >= 0 & km$incidence <= 1
    R <- rbind(R, km[i,,drop=FALSE])
  }
  if(onlydata) return(R)
  g <- ggplot(R, aes(x=.data$surv, y=1 - .data$incidence, color=.data$Statistic)) + geom_line() +
    geom_abline(intercept=0, slope=1, alpha=0.3) +
    xlab(xlab) + ylab(ylab) + guides(color=guide_legend(title=expression(y)))
  if(hare) g <- g + facet_wrap(~ .data$Type, nrow=nrow)
  g
}
