#' ordESS
#'
#' Ordinal Model Effective Sample Size
#'
#' For a standard ordinal model fitted with `orm`, returns the effective sample size (ESS) component of the `stats` part of the fit object if there were no censored data.  Otherwise `ordESS` assumes that `y=TRUE` and `lpe=TRUE` were given to `orm`, and an analysis of the effective sample size per censored observation is given, as a function of the censoring time, or in the case of interval censored data, o function of the width of the interval.
#'
#' @param fit a model fitted by `orm` with `y=TRUE, lpe=TRUE`
#'
#' @returns a `ggplot2` object
#' @md
#' @author Frank Harrell
#' @export
#'
#' @examples
#' \dontrun{
#' f <- orm(Ocens(y1, y2) ~ x, y=TRUE, lpe=TRUE)
#' ordESS(f)
#' }
ordESS <- function(fit) {
  if('y' %nin% names(fit) || 'lpe' %nin% names(fit) || ! length(fit$Ncens1) ||
     sum(fit$Ncens1) == 0) {
    message('Fit did not specify y=TRUE, lpe=TRUE.  Returning ESS stored with fit.')
    return(fit$stats['ESS'])
  }
  Nc    <- fit[['Ncens1']]
  Y     <- fit[['y']] 
  p     <- fit[['lpe']]
  v     <- fit$yunique
  ll    <- -2 * log(p)
  n     <- fit$stats['Obs']
  ncens <- sum(Nc)
  nun   <- n - ncens
  y1 <- as.vector(Y[, 1])
  y2 <- as.vector(Y[, 2])
  # Compute multiplier that makes sum of -2 LL for uncensored observations add up to
  # the number
  k     <- nun / sum(ll[is.finite(y1) & is.finite(y2) & (y1 == y2)])
  ll    <- ll * k
  # Use this multiplier to get per-observation effective sample size for each
  # uncensored observation
  d <- NULL
  g <- function(type)
        data.frame(dur = dur, ESS = ll[i],
                   type = paste0(type, ' (ESS=', round(sum(ll[i]), 1), ')'))
  ni <- Nc['interval']

  if(Nc['left'] > 0) {
    i <- is.infinite(y1)
    dur <- v[y2[i]]
    d <- rbind(d, g('Left Censored'))
  }
  if(Nc['right'] > 0) {
    i <- is.infinite(y2)
    dur <- y1[i]
    d <- rbind(d, g('Right Censored'))
  }
  if(ni > 0) {
    i <- is.finite(y1) & is.finite(y2) & (y1 < y2)
    dur <- (y2 - y1)[i]
    d <- rbind(d, g('Interval Censored'))
  }
  xl  <- if(ni > 0) 'Duration' else 'Censoring Point'
  cap <- if(ni > 0) 'Duration is censoring point or width of interval'
  ggplot(d) + aes(x=dur, y=ESS, color=type) + geom_point() + geom_smooth() +
         xlab(xl) + ylab('ESS Per Observation') +
         labs(caption=cap) +
         guides(color=guide_legend(title=''))
}
utils::globalVariables(c('ESS', 'type'))
