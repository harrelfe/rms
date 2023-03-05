npsurv <- function(formula, data=environment(formula), subset,
                   weights, na.action=na.delete, ...) {

  callenv <- parent.frame()
  w <- list(formula=formula, data=data, na.action=na.action)
  if(! missing(weights)) w$weights <- eval(substitute(weights), data, callenv)
  if(! missing(subset )) w$subset  <- eval(substitute(subset),  data, callenv)

  g <- do.call('model.frame', w)
  f <- do.call('survfit', w)

  f$maxtime     <- max(f$time)
  Y             <- g[[1]]
  f$units       <- units(Y)
  f$time.label  <- label(Y, type='time')
  f$event.label <- label(Y, type='event')

  strat <- rep('', NROW(Y))
  if(length(f$strata)) {
    X <- g[-1]
    nx <- ncol(X)
    for(j in 1 : nx)
      strat <- paste(strat, names(X)[j], '=', as.character(X[[j]]),
                     if(j < nx) ', ', sep='')
  }

  f$numevents <- if(inherits(f, 'survfitms')) {
    ## competing risk data; survfit.formula forgot to compute
    ## number of events for each state
    states <- attr(Y, 'states')
    state  <- factor(Y[, 'status'], 0 : length(states),
                     attr(Y, 'inputAttributes')$event$levels)                   #                                    c('censor', states))
    table(strat, state)
  }
  else tapply(Y[, 'status'], strat, sum, na.rm=TRUE)

  ## Compute person-time of exposure while we're at it
  f$exposure  <- tapply(Y[, 1], strat, sum, na.rm=TRUE)
  
  f$call <- match.call()
  class(f) <- c('npsurv', class(f))
  f
}
