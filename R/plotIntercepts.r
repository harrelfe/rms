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
#' @param dots set to `TRUE` to show solid dots at the intecept values
#' @param logt set to `TRUE` to use a log scale for the x-axis
#'
#' @returns `ggplot2` object
#' @export
#' @md
#' @author Frank Harrell
#'
#' @examples
#' \dontrun{
#' f <- orm(y ~ x1 + x2 + x3)
#' plotIntercepts(f)
#' }
plotIntercepts <- function(fit, dots=FALSE, logt=FALSE) {
  if(! inherits(fit, c('lrm', 'orm'))) stop('fit must be from lrm or orm')
  isorm <- inherits(fit, 'orm')

  # opar <- par(mar=c(4,4,2,3), mgp=c(3-.75,1-.5,0))
  # on.exit(par(opar))

  ns     <- num.intercepts(fit)
  alpha  <- coef(fit)[1 : ns]
  y      <- fit$yunique[-1]
  ylabel <- fit$ylabel
  yname  <- all.vars(fit$sformula)[1]
  if(! length(ylabel) || ylabel == '') ylabel <- yname
  if(isorm) ylabel <- fit$yplabel
  
  # plot(y, alpha, log=if(logt) 'x' else '',
  #      xlab=ylabel, ylab='Intercept', pch=20, cex=if(dots) 0.7 else 0)
  # segments(y[-ns], alpha[-ns], y[-1], alpha[-ns])                # horizontals
  # segments(y[-1],  alpha[-ns], y[-1], alpha[-1], col='gray85')   # verticals

  xtrans <- if(logt) 'log' else 'identity'
  npretty <- 10
  if(xtrans == 'identity') {
    xbreaks <- pretty(y, npretty)
    labels  <- format(xbreaks)
  } else if(FALSE) {
    xbreaks <- pretty(y, 2 * npretty)
    xbreaks <- xbreaks[xbreaks > 0]
    if(xbreaks[1] >= 1) xbreaks <- c(0.1, 0.25, 0.5, 0.75, xbreaks)
    if(xbreaks[1] > 0.5) xbreaks <- c(0.1, 0.25, 0.5, xbreaks)
    if(xbreaks[1] > 0.1) xbreaks <- c(0.1, xbreaks)
    lxb     <- log(xbreaks)
    lxbr    <- rmClose(lxb, 0.06)
    xbreaks <- xbreaks[lxb %in% lxbr]
  }

  g <- ggplot(mapping=aes(x=y, y=alpha)) + geom_step() +
         xlab(ylabel) + ylab(expression(alpha))
  if(logt) g <- g + scale_x_log10(guide='axis_logticks')
  else     g <- g + scale_x_continuous(breaks=xbreaks)

  if(dots) g <- g + geom_point(size=0.5)
  g
}
