#' Title Plot npsurv Nonparametric Survival Curves Using ggplot2
#'
#' @param data the result of npsurv
#' @param mapping unused
#' @param conf set to `"none"` to suppress confidence bands
#' @param trans the name of a transformation for the survival probabilities to use in drawing the y-axis scale.  The default is no transformation, and other choices are `"logit", "probit", "loglog"`.  `"loglog"` represents \eqn{-log(-log(S(t)))}
#' @param logt set to `TRUE` to use a log scale for the x-axis
#' @param curtail set to a (lower, upper) 2-vector to curtail survival probabilities and confidence limits before transforming and plotting
#' @param xlab x-axis label, the default coming from `fit`
#' @param ylab y-axis label, the default coming from `fit`
#' @param abbrev.label set to `TRUE` to abbreviate strata levels
#' @param levels.only set to `FALSE` to keep the original strata name in the levels
#' @param alpha transparency for confidence bands
#' @param facet when strata are present, set to `TRUE` to facet them rather than using colors on one panel
#' @param npretty the number of major tick mark labels to be constructed by [scales::breaks_pretty()] or [pretty()].  For transformed scales, twice this number is used.
#' @param onlydata set to `TRUE` to return the data frame to be plotted, and no plot
#' @param ... ignored
#' @param environment unused
#' @md
#' @author Frank Harrell
#' @returns a `ggplot2` object, if `onlydata=FALSE`
#' @method ggplot npsurv
#' @export
#'
#' @examples
#' set.seed(1)
#' g <- c(rep('a', 500), rep('b', 500))
#' y <- exp(-1 + 2 * (g == 'b') + rlogis(1000) / 3)
#' f <- npsurv(Surv(y) ~ g)
#' ggplot(f, trans='logit', logt=TRUE)

ggplot.npsurv <-
  function(data, mapping, conf=c('bands', 'none'),
           trans=c('identity', 'logit', 'probit', 'loglog'),
           logt=FALSE, curtail=c(0,1),
           xlab, ylab='Survival Probability',
           abbrev.label=FALSE, levels.only=TRUE,
           alpha=0.15, facet=FALSE, npretty=10, onlydata=FALSE, ..., environment) {

  fit      <- data
  trans    <- match.arg(trans)
  conf     <- match.arg(conf)
  conf     <- conf == 'bands' && length(fit$lower)

  units <- Punits(fit$units)

  if(missing(xlab))
    xlab <- if(length(fit$time.label) && fit$time.label != '')
              labelPlotmath(fit$time.label, units)
            else if(units != '') upFirst(units)
            else 'Time'

  slev  <- names(fit$strata)
  if(levels.only) slev <- gsub('.*=', '', slev)
  sleva <- if(abbrev.label) abbreviate(slev) else slev
  ns    <- length(slev)

  ap    <- trans == 'identity' && ! logt
  atime <- if(ap) 0  else numeric()
  asurv <- if(ap) 1  else numeric()
  alim  <- if(ap) NA else numeric()

  if(ns <= 1) {
    d <- data.frame(time = c(atime, fit$time),
                    surv = c(asurv, fit$surv))
    if(length(fit$lower)) {
      d$lower <- c(alim, fit$lower)
      d$upper <- c(alim, fit$upper)
    }
  }
  else {
    gr <- rep(slev, fit$strata)
     d <- NULL
     for(i in 1 : ns) {
      j <- gr == slev[i]
      dat <- data.frame(gr   = sleva[i],
                        time = c(atime, fit$time[j]),
                        surv = c(asurv, fit$surv[j]) )
      if(length(fit$lower)) {
        dat$lower <- c(alim, fit$lower[j])
        dat$upper <- c(alim, fit$upper[j])
      }
      d <- rbind(d, dat)
    }
    d$gr <- factor(d$gr, sleva)
  }

  if(trans == 'loglog')
    loglog <- function() trans_new('loglog', function(x) - log(-log(x)),
                                             function(x) exp(-exp(-x)),
                                           breaks = breaks_pretty(n=npretty))
  xtrans <- if(logt) 'log' else 'identity'
  if(! missing(curtail)) {
    curt <- function(x) pmin(curtail[2], pmax(x, curtail[1]))
    d$surv <- curt(d$surv)
    if(conf) {
      d$lower <- curt(d$lower)
      d$upper <- curt(d$upper)
    }
  }
  if(trans != 'identity') {
    d <- d[! is.na(d$surv) & d$surv > 0e0 & d$surv < 1e0, ]
    if(conf) {
      i <- (! is.na(d$lower) & d$lower == 0e0) | (! is.na(d$upper) & d$upper == 1e0)
      d$lower[i] <- NA
      d$upper[i] <- NA
    }
  }

  if(onlydata) return(d)

  pb <- breaks_pretty(n = npretty)
  if(xtrans == 'identity') {
    xbreaks <- pretty(d$time, npretty)
    labels  <- format(xbreaks)
  } else {
    xbreaks <- pretty(d$time, 2 * npretty)
    xbreaks <- xbreaks[xbreaks > 0]
    if(xbreaks[1] >= 1) xbreaks <- c(0.1, 0.25, 0.5, 0.75, xbreaks)
    lxb     <- log(xbreaks)
    lxbr    <- rmClose(lxb, 0.06)
    xbreaks <- xbreaks[lxb %in% lxbr]
  }

  if(trans == 'identity') ybreaks <- breaks_pretty(n = npretty)
  else {
    ybreaks <- pretty(d$surv, 2 * npretty)
    ybreaks <- sort(unique(c(ybreaks, seq(0.9, 0.99, by=0.01))))
    ybreaks <- ybreaks[ybreaks > 0 & ybreaks < 1]
    tyb <- switch(trans,
                  logit  = qlogis(ybreaks),
                  probit = qnorm(ybreaks),
                  loglog = -log(-log(ybreaks)) )
    tybr <- rmClose(tyb, 0.04)
    ybreaks <- ybreaks[tyb %in% tybr]
  }

  w <- list(if(ns > 1 && ! facet)
              geom_step(aes(color=.data$gr)) else geom_step(),
            if(conf && ns > 1 && ! facet)
              geom_stepconfint(aes(ymin=.data$lower, ymax=.data$upper, fill=.data$gr), alpha=alpha),
            if(conf && (facet || ns < 2))
              geom_stepconfint(aes(ymin=.data$lower, ymax=.data$upper), alpha=alpha),
            scale_x_continuous(transform=xtrans, breaks=xbreaks),  # pb),
            if(trans == 'identity')
              scale_y_continuous(breaks=ybreaks),
            if(trans %in% c('logit', 'probit'))
              scale_y_continuous(transform=trans, breaks=ybreaks),
            if(trans == 'loglog')
              scale_y_continuous(transform=loglog(), breaks=ybreaks),
            if(facet && ns > 1)
              facet_wrap(~ .data$gr),
            if(! facet && ns > 1)
              guides(color=guide_legend(title=''), fill=guide_legend(title=''))

  )
  ggplot(d, aes(x=.data$time, y=.data$surv)) + xlab(xlab) + ylab(ylab) + w
}
