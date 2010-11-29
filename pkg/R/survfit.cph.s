survfit.cph <- function(formula, newdata, se.fit=TRUE, conf.int=.95, 
                        individual=FALSE, type, vartype,
                        conf.type=c('log', 'log-log', 'plain', 'none'),
                        ...)
{
  require(survival)
  conf.type <- match.arg(conf.type)
  object <- formula
  
  if(!length(object$x)) stop('must use x=TRUE with fit')
  y <- object$y
  if(!length(y)) stop('must use y=TRUE with fit')

  strata <- attr(y, 'strata')
  Terms <- terms(formula)
  attr(Terms, 'specials')$strata <- attr(Terms, 'specials')$strat

#  object$strata <- strata
  smo <- NULL    # survival package model object for survfit.coxph
  object$terms  <- Terms
  object$n      <- sum(object$n)
  class(object) <- 'coxph'
  
  g <- list(formula=object)
  if(!missing(type   )) g$type    <- type
  if(!missing(vartype)) g$vartype <- vartype
  g$conf.type <- conf.type

  rq <- NULL
  if(!missing(newdata))
    {
      newdata <- predictrms(object, newdata, type='x', expand.na=FALSE)
      rq <- attr(newdata, 'strata')
      g$newdata <- newdata
    }
  survfits <- survfit
  g$censor <- FALSE  # don't output censored values

  g <- do.call('survfits', g)
  g$requested.strata <- rq
  class(g) <- c('survfit.cph', class(g))
  g$call <- object$call
  g
}

