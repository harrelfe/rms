calibrate.orm <- function(fit, 
                          method="boot", u, m=150, pred, B=40,
                          bw=FALSE, rule="aic",
                          type="residual", sls=.05, aics=0, force=NULL,
                          estimates=TRUE, pr=FALSE, what="observed-predicted",
                          val.surv.args=list(method='smoothkm', eps=30),
                          ...)
{
  call    <- match.call()

  deb <- Fdebug('calibrate.debug')

  unit <- fit$units
  if(unit=="") unit <- "unit"

  y       <- fit[['y']]
  N       <- NROW(y)
  Nevents <- N

  Ncens1 <- fit$Ncens1
  if(length(Ncens1)) {
    if(Ncens1['left'] + Ncens1['interval'] > 0)
      stop('left or interval censoring not implemented') 
    Nevents <- N - Ncens1['right']
    if(Nevents < 5) stop('must have >= 5 uncensored observations')
  }

  xb       <- fit$linear.predictors
  survival <- survest(fit, times=u, conf.int=0)$surv
  deb(survival)

  if(missing(pred)) {
    lim  <- datadist(survival)$limits[c('Low:prediction','High:prediction'),]
    pred <- seq(lim[1], lim[2], length=100)
    deb(pred)
  }
  
  distance <- function(x, y, fit, iter, u, fit.orig, what="observed",
                       pred, ...) {
    deb(c(iter,
          num.intercepts(fit),
          length(fit$Design),
          if(NCOL(y) == 2) sum(y[, 1] == y[, 2])));   deb(pred);  deb(x)

    # x is X*beta, y is an Ocens object or a plain vector
    # x uses first intercept.  Adjust here to use interceptRef which
    # works better for survest.

    # Don't compute accuracy when < 5 uncensored observations
    if(NCOL(y) == 2 && sum(y[, 1] == y[, 2]) < 5) return(NA)

    kint  <- fit$interceptRef
    alpha <- coef(fit)   # don't care if betas are at the end
    x     <- x - alpha[1] + alpha[kint]

    psurv <- survest(fit, linear.predictors=x, times=u, conf.int=0)$surv
    deb(psurv)

    pred.obs <- do.call(val.surv,
                        c(list(fit, S=y, u=u, est.surv=psurv,
                               pred=pred), val.surv.args))
    deb(unclass(pred.obs))
    dist <- if(what=='observed') pred.obs$actualseq
            else                 pred.obs$actualseq - pred
    
    if(iter == 0 && pr) print(pred.obs)
    if(iter == 0) structure(dist, keepinfo=list(pred.obs=pred.obs)) else
    dist
  }

  ofit <- quickRefit(fit, what='fitter', storevals=FALSE, compstats=FALSE)

  reliability <-
    predab.resample(fit, method=method,
                    fit=ofit, measure=distance,
                    pr=pr, B=B, bw=bw, rule=rule, type=type,  
                    u=u, m=m, what=what, sls=sls, aics=aics,
                    force=force, estimates=estimates,
                    pred=pred, allow.varying.intercepts=TRUE, ...)
  kept     <- attr(reliability, 'kept')
  keepinfo <- attr(reliability, 'keepinfo')
  n    <- reliability[, "n"]
  rel  <- reliability[, "index.corrected"]
  opt  <- reliability[, "optimism"]
  
  rel <- cbind(mean.optimism=opt, mean.corrected=rel, n=n)
  
  pred.obs              <- keepinfo$pred.obs
  calibrated            <- pred.obs$actualseq
  calibrated.corrected  <- calibrated - opt
    
  structure(cbind(pred=pred,
                  reliability[, c("index.orig", "training", "test"), drop=FALSE],
                  rel, calibrated=calibrated,
                  calibrated.corrected=calibrated.corrected  ),
                  predicted=survival, kept=kept,
                  class="calibrate", u=u, units=unit, n=N, d=Nevents,
                  p=length(fit$coefficients), m=m, B=B, what=what, call=call)
}
