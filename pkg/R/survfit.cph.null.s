if(FALSE)
  survfit.cph.null <-
  function(formula, newdata, se.fit=TRUE, conf.int=.95, individual=FALSE,
           type, vartype,
           conf.type=c('log', 'log-log', 'plain', 'none'), ...)
{
  require(survival)
  conf.type <- match.arg(conf.type)
  f <- formula
  y <- f$y
  if(!length(y)) stop("must use y=TRUE with fit")
  n <- nrow(y)

  Strata <- attr(y, 'strata')
  if(length(Strata))
    {
      n.all <- table(Strata) # patch bug in survfit.coxph.null
      storeTemp(n.all)       # survfit.coxph.null can't find here (why?)
    }
  Terms <- terms(formula)
  attr(Terms,'specials')$strata <- attr(Terms,'specials')$strat
  
  g <- list(formula=
            list(call=f$call, loglik=f$loglik, residuals=f$residuals,
                 n=sum(f$n), linear.predictors=rep(0, n),
                 method=f$method, assign=f$assign,
                 y=y, strata=Strata, terms=Terms, formula=f$formula))
  class(g$formula) <- 'coxph.null'

  if(!missing(type))      g$type      <- type
  if(!missing(vartype))   g$vartype   <- vartype
  g$conf.type <- conf.type
  
  survfits <- survfit
  g <- do.call('survfits', g)
  g$call <- f$call
  class(g) <- c('survfit.cph.null', class(g))
  g
}
