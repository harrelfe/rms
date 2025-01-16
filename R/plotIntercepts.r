#' Plot Intercepts
#'
#' Plots the step function corresponding to the intercepts in a `orm` or `lrm` model.  This can be thought
#' of as the link function of the covariate-adjusted empirical cumulative distribution function
#' (actually 1 - ECDF).  It is
#' also related to q-q plots.  For example, if a probit link function is an appropriate choice, and the
#' residuals actually had a normal distribution (not needed by the semiparametric ordinal model), the step
#' function of the intercepts would form a straight line.
#'
#' @param fit an `orm` or `lrm` fit object, usually with a numeric dependent variable having many levels
#'
#' @returns nothing; only plots
#' @export
#' @md
#' @author Frank Harrell
#'
#' @examples
#' \dontrun{
#' f <- orm(y ~ x1 + x2 + x3)
#' plotIntercepts(f)
#' }
plotIntercepts <- function(fit) {
  if(! inherits(fit, 'lrm') && ! inherits(fit, 'orm')) stop('fit must be from lrm or orm')

  opar <- par(mar=c(4,4,2,3), mgp=c(3-.75,1-.5,0))
  on.exit(par(opar))

  ns     <- num.intercepts(fit)
  alpha  <- coef(fit)[1 : ns]
  y      <- fit$yunique[-1]
  yname  <- all.vars(fit$sformula)[1]

  plot(y, alpha, xlab=yname, ylab='Intercept')
  segments(y[-ns], alpha[-ns], y[-1], alpha[-ns])                # horizontals
  segments(y[-1],  alpha[-ns], y[-1], alpha[-1], col='gray85')   # verticals
}
