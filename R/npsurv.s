npsurv <- function(formula, data, subset, na.action, ...)
{
  m <- match.call(expand.dots = FALSE)
  m$... <- NULL
  m[[1]] <- as.name('survfit')
  m$formula <- formula
  f <- eval(m, sys.parent())
  m[[1]] <- as.name('model.frame')
  g <- eval(m, sys.parent())
  
  f$maxtime <- max(f$time)
  
  Y <- model.extract(g, 'response')
  f$units <- units(Y)
  f$time.label  <- label(Y, type='time')
  f$event.label <- label(Y, type='event')
  f$call <- match.call()
  class(f) <- c('npsurv', 'survfit')
  f
}
