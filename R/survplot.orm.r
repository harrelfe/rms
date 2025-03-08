#' Title Survival Curve Plotting
#'
#' Plots predicted survival curves with easy specification of predictor settings, with optional confidence bands.  For `orm` fits these are step functions, and for `psm` fits they are smooth curves.
#'
#' @param fit a fit produced by [orm()]; also works for [psm()] fits
#' @param ... list of factors with names used in model. The first factor listed is the factor used to determine different survival curves.  Any other factors are used to specify single constants to be adjusted to, when defaults given to fitting routine (through `limits`) are not used.  The value given to factors is the original coding of data given to fit, except that for categorical factors the text string levels may be specified.  The form of values given to the first factor are none (omit the equal sign to use default range or list of all values if variable is discrete), `"text"` if factor is categorical, `c(value1, value2, \dots)`, or a function which returns a vector, such as `seq(low,high,by=increment)`.  Only the first factor may have the values omitted.  In this case the `Low effect`, `Adjust to`, and `High effect` values will be used from `datadist` if the variable is continuous.  For variables not defined to `datadist`, you must specify non-missing constant settings (or a vector of settings for the one displayed variable).
#' @param xlab character string label for x-axis; uses the `plotmath`-style `yplabel` for the `y` variable stored in the fit if `xlab` is absent
#' @param ylab y-axis label, defaulting to `"Survival Probability"`
#' @param conf.int defaults to `FALSE` (same as specifying `0`); specify a positive value less than 1 to get two-sided confidence intervals utilizing approximate normality of linear predictors
#' @param conf not currently used
#' @param facet set to `TRUE` to have the first varying variable appear as a facet instead of as different colored step functions
#' @param nrow when faceting on one varying variable using `facet_wrap` specifies the number of rows to create
#' @param alpha transparency for confidence bands
#' @param adj.subtitle set to `FALSE` to not show a caption with the values of non-varying values (adjustment variables)
#' @param onlydata set to `TRUE` to return the data used in `ggplot2` plotting instead of the graphics object
#'
#' @returns if `onlydata` is left at its default value, a `ggplot2` graphics object for which additional layers may later be added
#' @seealso [Hmisc::geom_stepconfint()]
#' @export
#' @author Frank Harrell
#' md
#' @examples
#' set.seed(1)
#' d <- expand.grid(x1=c('a', 'b', 'c'), x2=c('A','B'), x3=1:2, irep=1:20)
#' y <- sample(1:10, nrow(d), TRUE)
#' dd <- datadist(d); options(datadist='dd')
#' f <- orm(y ~ x1 + x2 + x3, data=d)
#'
#' survplot(f, x1='a')
#' survplot(f, x1='a', conf.int=.95)
#' survplot(f, x1=c('a','b'), x2='A')
#' survplot(f, x1=c('a', 'b'), x2='A', conf.int=.95)
#' survplot(f, x1=c('a','b'), x2='A', facet=TRUE)
#' survplot(f, x1=c('a','b'), x2='A', facet=TRUE, conf.int=.95)
#'
#' survplot(f, x1=c('a', 'b'), x2=c('A', 'B'))
#' survplot(f, x1=c('a', 'b'), x2=c('A', 'B'), conf.int=.95)
#' survplot(f, x1=c('a', 'b'), x2=c('A', 'B'), facet=TRUE)
#'
#' survplot(f, x1=c('a', 'b'), x2=c('A', 'B'), x3=1:2)
#'
#' g <- psm(Surv(y) ~ x1 + x2 + x3, data=d)
#' survplot(g, x1=c('a','b'), x2=c('A', 'B'), ggplot=TRUE)  # calls survplot.orm
#' # See https://hbiostat.org/rmsc/parsurv#sec-parsurv-assess
#' # where nonparametric and parametric estimates are combined into one ggplot
#' options(datadist=NULL)

survplot.orm <-
  function(fit, ...,
           xlab, ylab='Survival Probability',
           conf.int=FALSE, conf=c("bands", "bars"),
           facet=FALSE, nrow=NULL, alpha=0.15,
           adj.subtitle=TRUE, onlydata=FALSE)
{
  conf      <- match.arg(conf)
  if(conf == 'bars') stop('conf="bars" is not yet implemented')
  ispsm     <- inherits(fit, 'psm')

  db        <- getOption('rmsdebug', FALSE)

  if(missing(xlab)) xlab <- fit$yplabel
  if(missing(ylab)) ylab <- 'Survival Probability'

  if(is.logical(conf.int)) {
    if(conf.int) conf.int <- .95	else conf.int <- 0
  }

  xadj    <- Predict(fit, type='model.frame', np=5,
                     factors=rmsArgs(substitute(list(...))))

  info    <- attr(xadj, 'info')
  varying <- info$varying
  nv      <- length(varying)
  nf      <- if(facet) nv else nv - 1  # no. of facet variables
  if(nf > 2)
    stop('cannot facet more than 2 varying predictors')
  use.color <- (nv > 0) && ! facet

  adjust <- if(adj.subtitle) info$adjust
  if(db) prn(llist(varying, nv, nf, use.color, xadj))
  nc     <- if(length(xadj)) nrow(xadj) else 1

  d <- NULL
  for(i in 1 : nc) {
    adj <- xadj[i, , drop=FALSE]
    w   <- survest(fit, newdata=adj, conf.int=conf.int)
    if(ispsm)
      w <- with(w, if(conf.int) data.frame(time, surv, lower, upper)
                           else data.frame(time, surv))
    # cbind complains about short row names from adj[varying]
    d   <- if(nv) suppressWarnings(rbind(d, cbind(adj[varying], w))) else rbind(d, w)
  }
  if(onlydata) return(d)

if(nv > 0) {
  cname <- names(d)[1]
  vs <- paste0('v', 1 : nv)
  names(d)[1 : nv] <- if(use.color) c('v', if(nv > 1) vs[1 : (nv - 1)]) else vs
  if(use.color) d$v <- factor(d$v)   # for ggplot2 color
}
g <- ggplot(d, aes(x=.data$time, y=.data$surv))
w <- list(
      if(use.color)   if(ispsm) geom_line(aes(color=.data$v)) else geom_step(aes(color=.data$v)),
      if(! use.color) if(ispsm) geom_line() else geom_step(),
      if(conf.int &&   use.color)
        if(ispsm) geom_ribbon(aes(ymin=.data$lower, ymax=.data$upper, fill=.data$v), alpha=alpha)
        else      geom_stepconfint(aes(ymin=.data$lower, ymax=.data$upper, fill=.data$v), alpha=alpha),
      if(conf.int && ! use.color)
        if(ispsm) geom_ribbon(aes(ymin=.data$lower, ymax=.data$upper), alpha=alpha)
        else      geom_stepconfint(aes(ymin=.data$lower, ymax=.data$upper), alpha=alpha),
      if(nf == 1) facet_wrap(~ .data$v1),
      if(nf == 2) facet_grid(.data$v1 ~ .data$v2),
      labs(x = xlab, y=ylab),
      if(use.color) guides(color=guide_legend(title=cname), fill=guide_legend(title=cname)),
      if(length(adjust)) labs(caption=paste('Adjusted to', as.character(adjust))),
      scale_x_continuous(expand = c(0.02, 0)) )
g + w
}
