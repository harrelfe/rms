#' Check Parallelism Assumption of Ordinal Semiparametric Models
#'
#' `orm` models are refitted as a series of binary models for a sequence of cutoffs
#' on the dependent variable.  Regression coefficients from this sequence are plotted
#' against cutoffs using `ggplot2` with one panel per regression coefficient.
#' When censoring is present, whether or not Y is
#' greater than or equal to the current cutoff is not always possible, and such
#' observations are ignored.
#' 
#' Whenver a cut gives rise to extremely high standard error for a regression coefficient,
#' the confidence limits are set to `NA`.  Unreasonable standard errors are determined from
#' the confidence interval width exceeding 7 times the standard error at the middle Y cut.
#' 
#' @param fit a fit object from `orm` with `x=TRUE, y=TRUE` in effect
#' @param which specifies which columns of the design matrix are assessed.  By default, all columns are analyzed.
#' @param terms set to `TRUE` to collapse all components of each predictor into a single column weighted by the original regression coefficients but scaled according to `scale`.  This means that each predictor will have a regression coefficient of 1.0 when refitting the original model on this transformed X matrix, before any further scaling.  Plots will then show the relative effects over time, i.e., the slope of these combined columns over cuts on Y, so that deviations indicate non-parallelism.  But since in this case only relative effects are shown, a weak predictor may be interpreted as having an exagerrated y-dependency if `scale='none'`.  `terms` detauls to `TRUE` when `onlydata=TRUE`.
#' @param m the lowest cutoff is chosen as the first Y value having at meast `m` observations to its left, and the highest cutoff is chosen so that there are at least `m` observations tot he right of it.  Cutoffs are equally spaced between these values.  If omitted, `m` is set to the minimum of 50 and one quarter of the sample size.
#' @param maxcuts the maximum number of cutoffs analyzed
#' @param lp plot the effect of the linear predictor across cutpoints instead of analyzing individual predictors
#' @param onlydata set to `TRUE` to return a data frame suitable for modeling effects of cuts, instead of constructing a graph.  The returned data frame has variables `Ycut, Yge_cut, obs`, and the original names of the predictors.  `Ycut` has the cutpoint on the original scale.  `Yge_cut` is `TRUE/FALSE` dependent on whether the Y variable is greater than or equal to `Ycut`, with `NA` if censoring prevented this determination.  The `obs` variable is useful for passing as the `cluster` argument to [robcov()] to account for the high correlations in regression coefficients across cuts.  See the example which computes Wald tests for parallelism where the `Ycut` dependence involves a spline function.  But since `terms` was used, each predictor is reduced to a single degree of freedom.
#' @param scale applies to `terms=TRUE`; set to `'none'` to leave the predictor terms scaled by regression coefficient so the coefficient of each term in the overall fit is 1.0.  The default is to scale terms by the interquartile-range (Gini's mean difference if IQR is zero) of the term.  This prevents changes in weak predictors over different cutoffs from being impressive.
#' @param conf.int confidence level for computing Wald confidence intervals for regression coefficients.  Set to 0 to suppress confidence bands.
#' @param alpha saturation for confidence bands
#'
#' @returns `ggplot2` object or a data frame
#' @export
#' @md
#' @author Frank Harrell
#' @examples
#' \dontrun{
#' f <- orm(..., x=TRUE, y=TRUE)
#' ordParallel(f, which=1:5)  # first 5 betas
#' 
#' getHdata(nhgh)
#' set.seed(1)
#' nhgh$ran <- runif(nrow(nhgh))
#' f <- orm(gh ~ rcs(age, 4) + ran, data=nhgh, x=TRUE, y=TRUE)
#' ordParallel(f)  # one panel per parameter (multiple parameters per predictor)
#' dd <- datadist(nhgh); options(datadist='dd')
#' ordParallel(f, terms=TRUE)
#' d <- ordParallel(f, maxcuts=30, onlydata=TRUE)
#' dd2 <- datadist(d); options(datadist='dd2')  # needed for plotting
#' g <- orm(Yge_cut ~ (age + ran) * rcs(Ycut, 4), data=d, x=TRUE, y=TRUE)
#' h <- robcov(g, d$obs)
#' anova(h)
# # Plot inter-quartile-range (on linear predictor "terms") age
# # effect vs. cutoff y
#' qu <- quantile(d$age, c(1, 3)/4)
#' qu
#' cuts <- sort(unique(d$Ycut))
#' cuts
#' z <- contrast(h, list(age=qu[2], Ycut=cuts),
#'                  list(age=qu[1], Ycut=cuts))
#' z <- as.data.frame(z[.q(Ycut, Contrast, Lower, Upper)])
#' ggplot(z, aes(x=Ycut, y=Contrast)) + geom_line() +
#'   geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2)
#' }
ordParallel <- function(fit, which, terms=onlydata, m, maxcuts=75, 
                        lp=FALSE, onlydata=FALSE, 
                        scale=c('iqr', 'none'),
                        conf.int=0.95, alpha=0.15) {

  scale <- match.arg(scale)
  Y <- fit[['y']]
  if(! length(Y)) stop('requires y=TRUE specified to orm')
  X <- fit[['x']]
  if(! lp && ! length(X)) stop('requires x=TRUE specified to orm')
  n       <- NROW(Y)
  isocens <- NCOL(Y) == 2
  if(isocens) {
    Y    <- Ocens2ord(Y)
    vals <- attr(Y, 'levels')
    YO   <- extractCodedOcens(Y, what=4, ivalues=TRUE)
    Y    <- YO$y
  }
  else Y <- recode2integer(Y, ftable=FALSE)$y - 1
  k      <- num.intercepts(fit)
  fitter <- quickRefit(fit, what='fitter')
  cfreq  <- cumsum(fit$freq)
  # if(max(cfreq) != n) stop('program logic error in ordParallel ', max(cfreq), ' ', n)
  if(missing(m)) m <- min(ceiling(n / 4), 50)
  if(n < 2 * m) stop('must have at least 2*m observations for m=', m)
  # Find first intercept with at least m observations preceeding it, and the
  # last one with at least m observations after it
  if(! isocens) vals <- as.numeric(names(cfreq))
  lev  <- 0 : k
  low  <- min(lev[cfreq >= m])
  high <- max(lev[cfreq <= n - m])
  ks   <- unique(round(seq(low, high, length=maxcuts)))

  zcrit <- qnorm((conf.int + 1) / 2)
  ylab <- expression(beta)

  if(lp) {
    X <- fit$linear.predictors
    if(! length(X)) stop('fit does not have linear.predictors')
    ylab <- 'Relative Coefficient'
    R   <- NULL
    for(ct in ks) {
      y <- if(isocens) geqOcens(YO$a, YO$b, YO$ctype, ct) else Y >= ct
      j <- ! is.na(y)
      f <- fitter(X, y, subset=j)
      co <- coef(f)[-1]
      v <- vcov(f, intercepts='none')
      s <- sqrt(diag(v))[which]
      w <- data.frame(cut=ct, beta=co)
      if(conf.int > 0) {
        w$lower <- co - zcrit * s
        w$upper <- co + zcrit * s
        }
      R <- rbind(R, w)
    }
    pr  <- pretty(vals[ks + 1], n=20)
    i   <- approx(vals, 0 : k, xout=pr, rule=2)$y
    ig  <- rmClose(i, minfrac=0.085 / 2)
    pr  <- pr[i %in% ig]
    g <- ggplot(R, aes(x=.data$cut, y=.data$beta)) + geom_line(alpha=0.35) + geom_smooth() +
           scale_x_continuous(breaks=ig, labels=format(pr)) +
           xlab(fit$yplabel) + ylab(ylab) +
           labs(caption=paste(length(ks), 'cuts,', m, 'observations beyond outer cuts'),
                subtitle=paste(fit$family, 'family'))
    if(conf.int > 0) g <- g + geom_ribbon(aes(ymin=.data$lower, ymax=.data$upper), alpha=alpha)
    return(g)
  }

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
  O     <- data.frame(x=names(co), beta=co)  # overall estimates
  if(conf.int > 0) {
    v  <- vcov(fit, intercepts='none')
    s  <- sqrt(diag(v))[which]
    O$lower <- co - zcrit * s
    O$upper <- co + zcrit * s
  }
  R   <- NULL
  D   <- list()
  obs <- 1 : n
  ic  <- 0
  mid <- ks[which.min(abs(ks - median(ks)))]   # cut closest to middle cut

  for(ct in ks) {
    y <- if(isocens) geqOcens(YO$a, YO$b, YO$ctype, ct) else Y >= ct
    j <- ! is.na(y)
    if(onlydata) {
      ic <- ic + 1
      XX <- X[j,,drop=FALSE]
      D[[ic]] <- data.frame(Ycut=vals[ct + 1], XX, Yge_cut=y[j], obs=obs[j])
    } else {
      f <- if(terms) fitter(X, y, subset=j) else fitter(y=y, subset=j)
      co <- coef(f)[-1]
      co <- co[which]
      v <- vcov(f, intercepts='none')
      s <- sqrt(diag(v))[which]
      w <- data.frame(cut=ct, x=names(co), beta=co)
      if(conf.int > 0) {
        if(ct == mid) secm <- structure(s, names=names(co))
        w$lower <- co - zcrit * s
        w$upper <- co + zcrit * s
        }
      R <- rbind(R, w)
    }
  }

  # If any standard errors blew up, set confidence limits to NA for those
  # Test: confidence interval width > 7 times the standard error at the middle cut
  if(! onlydata && conf.int > 0) {
    for(xnam in names(co)) {
      i <- which(R$x == xnam)
      j <- R[i, 'upper'] - R[i, 'lower'] > 7 * secm[xnam]
      if(any(j)) {
        R[i[j], 'upper'] <- NA
        R[i[j], 'lower'] <- NA
      } 
    }
  }

  if(onlydata) {
    D <- do.call(rbind, D)
    return(D)
  }

  R$x <- factor(R$x, names(co), names(co))
  O$x <- factor(O$x, names(co), names(co))
  pr  <- pretty(vals[ks + 1], n=20)
  i   <- approx(vals, 0 : k, xout=pr, rule=2)$y
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
