# Use x= if input is a design matrix, newdata= if a data frame or data matrix
# or vector.

survest.orm <- function(fit, newdata, linear.predictors, x, times, fun, 
                        loglog=FALSE, conf.int=0.95,
                        what=c("survival", "parallel"), 
                        ...) { 

  what <- match.arg(what)
  if(what=='parallel') conf.int <- FALSE

  S <- Survival(fit)
  
  if(missing(fun)) fun <- if(loglog) function(x) logb(ifelse(x == 0 | x == 1, NA, x))
  else function(x) x

  if(conf.int > 0 && ! missing(linear.predictors)) {
    warning('conf.int set to 0 since linear.predictors specified')
    conf.int <- 0
  }

  p <- length(fit$coef) - num.intercepts(fit)
  
  if(missing(linear.predictors)) {
    linear.predictors <- if(p == 0) 0
    if(p > 0 && missing(x) && missing(newdata)) {
      if(missing(times)) stop('specify times= if using linear predictors from fit')
      linear.predictors <- fit$linear.predictors
      if(conf.int > 0)
        stop("may not specify conf.int unless x or newdata given")
      rnam <- names(linear.predictors)
    }
    else {
      if(missing(x)) x <- predict(fit, newdata, type="x")
      rnam <- dimnames(x)[[1]]
    }
  }
  else  rnam <- names(linear.predictors)
  
  if(what == 'parallel') {
    if(length(times) > 1 && (length(times) != length(linear.predictors)))
      stop('length of times must = 1 or number of subjects when what="parallel"')
    return(S(times, linear.predictors))
  }
  
  if(missing(times)) times <- NULL
  if(missing(x) && missing(newdata)) x <- NULL
  nt <- length(times)
  n <- length(linear.predictors)
  
  if(n > 1 & ! length(times))
    warning("should specify times if getting predictions for >1 obs.")
  
  # surv <- drop(outer(linear.predictors, times, FUN=comp, Trans=trans))

  S(times, linear.predictors, X=x, conf.int=conf.int)
  }