residuals.cph <-
  function(object,
           type = c("martingale", "deviance", "score",
             "schoenfeld", "dfbeta", "dfbetas", "scaledsch","partial"), ...)
  {
    type <- match.arg(type)
    x <- object[['x']]
    y <- object[['y']]
    if(type != 'martingale' && !length(x))
      stop('you must specify x=TRUE in the fit')
    if(type %nin% c('deviance','martingale') && !length(y))
      stop('you must specify y=TRUE in the fit')

    strata <- attr(x, 'strata')
    if(length(strata))
      {
      object$strata <- strata
      terms <- terms(object)
      attr(terms,'specials')$strata <- attr(terms,'specials')$strat
      object$terms <- terms
    }
    survival:::residuals.coxph(object, type=type, ...)
  }
