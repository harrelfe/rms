survest <- function(fit, ...) UseMethod("survest")

##Use x= if input is a design matrix, newdata= if a data frame or data matrix
#or vector.  Can specify (centered) linear predictor values instead (linear.predictors).
#Strata is attached to linear.predictors or x as "strata" attribute.
#data matrix assumes that categorical variables are coded with integer codes

survest.psm <- function(fit, newdata, linear.predictors, x, times, fun, 
                        loglog=FALSE, conf.int=.95,
                        what=c("survival","hazard","parallel"), 
                        ...) {   # ... so survplot will work

  what <- match.arg(what)
  if(what=='parallel') conf.int <- FALSE

  trans <- switch(what, survival=Survival(fit), hazard=Hazard(fit),
                  parallel=Survival(fit))

  if(missing(fun)) fun <- if(loglog) function(x) logb(ifelse(x==0|x==1,NA,x))
  else function(x) x

  if(what=="hazard" & conf.int>0) {
    warning('conf.int ignored for what="hazard"')
    conf.int <- FALSE
  }

  if(conf.int > 0) {
    cov <- vcov(fit, regcoef.only=TRUE)  # ignore scale
    if(!missing(linear.predictors))	{
      warning("conf.int set to 0 since linear.predictors specified")
      conf.int <- 0
    }
  }

  if(any(attr(fit,'class')=='pphsm'))
    stop("fit should not have been passed thru pphsm")

  nvar <- length(fit$coef) - num.intercepts(fit)
  
  if(missing(linear.predictors)) {
    if(nvar > 0 && missing(x) && missing(newdata)) {
      linear.predictors <- fit$linear.predictors
      if(conf.int > 0)
        stop("may not specify conf.int unless x or newdata given")
      rnam <- names(linear.predictors)
    }
    else {
      if(nvar==0) {
        x <- as.matrix(1)	# no predictors
        linear.predictors <- fit$coef[1]
      } else {
        if(missing(x)) x <- cbind(Intercept=1, predict(fit, newdata, type="x"))
        linear.predictors <- matxv(x, fit$coef)
      }
      if(conf.int > 0) {
        g1 <- drop(((x %*% cov) * x) %*% rep(1, ncol(x)))
        last <- {
          nscale <- length(fit$icoef) - 1
          ncol(fit$var) - (1 : nscale) + 1
        }
        g2 <- drop(x %*% fit$var[-last, last, drop=FALSE])
      }
      rnam <- dimnames(x)[[1]]
    }
  }
  else  rnam <- names(linear.predictors)
  
  if(what == 'parallel') {
    if(length(times)>1 && (length(times) != length(linear.predictors)))
      stop('length of times must = 1 or number of subjects when what="parallel"')
    return(trans(times, linear.predictors))
  }
  
  if(missing(times)) times <- seq(0, fit$maxtime, length=200)
  nt <- length(times)
  n <- length(linear.predictors)
  
  if(n > 1 & missing(times))
    warning("should specify times if getting predictions for >1 obs.")
  
  if(conf.int>0) zcrit <- qnorm((conf.int + 1) / 2)
  
  comp <- function(a, b, Trans) Trans(b, a)
  surv <- drop(outer(linear.predictors, times, FUN=comp, Trans=trans))
  
  if(conf.int > 0 && (nt==1 || n==1)) {
    dist <- fit$dist
    link <- survreg.distributions[[dist]]$trans
    z    <- if(length(link)) link(times) else times
    sc   <- fit$scale   ## TODO: generalize for vector of scale params
    logtxb <- outer(linear.predictors, z, function(a,b) b - a)
    se   <- sqrt(g1 + logtxb * (2 * g2 + logtxb * fit$var[last, last])) / sc
    prm  <- 0
    tm   <- if(length(link)) 1 else 0
    
    lower <- trans(tm,-drop(logtxb / sc + zcrit * se), parms=prm)
    upper <- trans(tm,-drop(logtxb / sc - zcrit * se), parms=prm)
    if(what=='survival') {
      lower[times == 0] <- 1
      upper[times == 0] <- 1
    }
    std.err <- drop(se)
  }
  
  if(nt==1 | n==1) {
    surv <- fun(surv); surv[is.infinite(surv)] <- NA
      if(conf.int > 0) {
        lower <- fun(lower); lower[is.infinite(lower)] <- NA
        upper <- fun(upper); upper[is.infinite(upper)] <- NA
        retlist <- list(time=times,surv=surv,
                        lower=lower,upper=upper,
                        std.err=std.err,
                        linear.predictors=linear.predictors)
      }
      else retlist <- list(time=times,surv=surv,
                           linear.predictors=linear.predictors)
    retlist <- structure(c(retlist,
                           list(conf.int=conf.int, units=fit$units,
                                  n.risk=fit$stats["Obs"],
                                  n.event=fit$stats["Events"], what=what)),
                           class='survest.psm')
      return(retlist)
    }
  
  if(n==1) names(surv) <- format(times) else {
    if(is.matrix(surv))
      dimnames(surv) <- list(rnam, format(times))
    else names(surv) <- rnam
  }
  
  surv
}


print.survest.psm <- function(x, ...)
{
  cat('\nN:',x$n.risk,'\tEvents:',x$n.event)
  z <- if(length(unique(x$time)) > 1) data.frame(Time=x$time) else {
    cat('\tTime:',x$time[1],' ',x$units,'s',sep='')
    data.frame(LinearPredictor=x$linear.predictors)
  }
  cat('\n\n')
  z$whatever  <- x$surv
  names(z)[2] <- x$what
  if(x$conf.int)
    {
      z$Lower <- x$lower
      z$Upper <- x$upper
      z$SE    <- x$std.err
    }
  print.data.frame(z)
  invisible()
}
