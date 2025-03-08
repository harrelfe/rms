#' Check Parallelism Assumption of Ordinal Semiparametric Models
#'
#' `orm` models are refitted as a series of binary models for a sequence of cutoffs
#' on the dependent variable.  Regression coefficients from this sequence are plotted
#' against cutoffs using `ggplot2` with one panel per regression coefficient.
#' When censoring is present, whether or not Y is
#' greater than or equal to the current cutoff is not always possible, and such
#' observations are ignored.
#' @param fit a fit object from `orm` with `x=TRUE, y=TRUE` in effect
#' @param which specifies which columns of the design matrix are assessed.  By default, all columns are analyzed.
#' @param terms set to `TRUE` to collapse all components of each predictor into a single column weighted by the original regression coefficients but scaled according to `scale`.  This means that each predictor will have a regression coefficient of 1.0 when refitting the original model on this transformed X matrix, before any further scaling.  Plots will then show the relative effects over time, i.e., the slope of these combined columns over cuts on Y, so that deviations indicate non-parallelism.  But since in this case only relative effects are shown, a weak predictor may be interpreted as having an exagerrated y-dependency if `scale='none'`.
#' @param m the lowest cutoff is chosen as the first Y value having at meast `m` observations to its left, and the highest cutoff is chosen so that there are at least `m` observations tot he right of it.  Cutoffs are equally spaced between these values.  If omitted, `m` is set to the minimum of 50 and one quarter of the sample size.
#' @param maxcuts the maximum number of cutoffs analyzed
#' @param scale applies to `terms=TRUE`; set to `'none'` to leave the predictor terms scaled by regression coefficient so the coefficient of each term in the overall fit is 1.0.  The default is to scale terms by the interquartile-range (Gini's mean difference if IQR is zero) of the term.  This prevents changes in weak predictors over different cutoffs from being impressive.
#' @param conf.int confidence level for computing Wald confidence intervals for regression coefficients.  Set to 0 to suppress confidence bands.
#' @param alpha saturation for confidence bands
#'
#' @returns `ggplot2` object
#' @export
#' @md
#' @author Frank Harrell
#' @examples
#' \dontrun{
#' f <- orm(..., x=TRUE, y=TRUE)
#' ordParallel(f, which=1:5)  # first 5 betas
#' }
ordParallel <- function(fit, which, terms=FALSE, m, maxcuts=75, scale=c('iqr', 'none'),
                        conf.int=0.95, alpha=0.15) {

  scale <- match.arg(scale)
  Y <- fit[['y']]
  if(! length(Y)) stop('requires y=TRUE specified to orm')
  X <- fit[['x']]
  if(! length(X)) stop('requires x=TRUE specified to orm')
  isocens <- NCOL(Y) == 2
  if(isocens) {
    Y  <- Ocens2ord(Y)
    YO <- extractCodedOcens(Y, what=4, ivalues=TRUE)
    Y  <- YO$y
  }
  else Y <- recode2integer(Y, ftable=FALSE)$y - 1
  k <- num.intercepts(fit)
  fitter <- quickRefit(fit, what='fitter')
  cfreq  <- cumsum(fit$freq)
  n <- max(cfreq)
  if(missing(m)) m <- min(ceiling(n / 4), 50)
  if(n < 2 * m) stop('must have at least 2*m observations for m=', m)
  # Find first intercept with at least m observations preceeding it, and the
  # last one with at least m observations after it
  lev <- as.numeric(names(cfreq))
  lev  <- 0 : k
  low  <- min(lev[cfreq >= m])
  high <- max(lev[cfreq <= n - m])
  ks   <- unique(round(seq(low, high, length=maxcuts)))

  zcrit <- qnorm((conf.int + 1) / 2)
  ylab <- expression(beta)
  if(terms) {
    X    <- predict(fit, type='terms')
    ylab <- 'Relative Coefficient'
    if(scale == 'iqr') {
      ylab <- 'Effect in IQR Units'
      iqr <- function(x) {
        d <- diff(quantile(x, c(0.25, 0.75)))
        if(d == 0e0) d <- GiniMd(d)
        d
      }
      iq <- apply(X, 2, iqr)
      X <- sweep(X, 2, iq, '/')
    }
    fit <- fitter(X)
  }
  co    <- coef(fit)[-(1 : k)]
  p     <- length(co)
  if(missing(which)) which <- 1 : p
  co    <- co[which]
  O     <- data.frame(x=names(co), beta=co)
  if(conf.int > 0) {
    v  <- vcov(fit, intercepts='none')
    s  <- sqrt(diag(v))[which]
    O$lower <- co - zcrit * s
    O$upper <- co + zcrit * s
  }
  R <- NULL
  for(ct in ks) {
    y <- if(isocens) geqOcens(YO$a, YO$b, YO$ctype, ct) else Y >= ct
    j <- ! is.na(y)
    f <- if(terms) fitter(X, y, subset=j) else fitter(y=y, subset=j)
    co <- coef(f)[-1]
    co <- co[which]
    v <- vcov(f, intercepts='none')
    s <- sqrt(diag(v))[which]
    w <- data.frame(cut=ct, x=names(co), beta=co)
    if(conf.int > 0) {
      w$lower <- co - zcrit * s
      w$upper <- co + zcrit * s
    }
    R <- rbind(R, w)
  }
  R$x <- factor(R$x, names(co), names(co))
  O$x <- factor(O$x, names(co), names(co))
  lev <- fit$yunique
  pr  <- pretty(lev[ks + 1], n=20)
  i   <- approx(lev, 0 : k, xout=pr, rule=2)$y
  ig  <- rmClose(i, minfrac=0.085 / ifelse(length(which) == 1, 2, 1))
  pr  <- pr[i %in% ig]
  g <- ggplot(R, aes(x=.data$cut, y=.data$beta)) + geom_line(alpha=0.35) + geom_smooth() +
        facet_wrap(~ .data$x, scales=if(terms && scale=='iqr')'fixed' else 'free_y') +
        scale_x_continuous(breaks=ig, labels=format(pr)) +
        xlab(fit$yplabel) + ylab(ylab) +
        labs(caption=paste(length(ks), 'cuts,', m, 'observations beyond outer cuts'),
             subtitle=paste(fit$family, 'family'))
  if(conf.int > 0) g <- g + geom_ribbon(aes(ymin=.data$lower, ymax=.data$upper), alpha=alpha)
  g <- g + geom_hline(aes(yintercept = .data$beta), data=O, color='red', alpha=0.4)
  if(conf.int > 0)
    g <- g + geom_segment(aes(x=ct[1], y=.data$lower,
                          xend=ct[1], yend=.data$upper),
                          data=O, color='red', alpha=0.4)
  g
}
