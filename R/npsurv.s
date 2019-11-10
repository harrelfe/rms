npsurv <- function(formula, data, subset, na.action, ...)
{
  M <- match.call()
  m <- M
  m[[1]] <- as.name('model.frame')
  m[names(m) %nin% c('', 'formula', 'data', 'subset', 'na.action')] <- NULL
  g <- eval(m, sys.parent())
  Y <- model.extract(g, 'response')

  m <- M
  m[[1]] <- as.name('survfit')
  m$formula <- formula
  f <- eval(m, sys.parent())
  f$maxtime     <- max(f$time)
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
