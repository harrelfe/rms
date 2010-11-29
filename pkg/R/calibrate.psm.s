calibrate.psm <- function(fit, cmethod=c('hare', 'KM'),
                          method="boot", u, m=150, pred, cuts, B=40,
                          bw=FALSE, rule="aic",
                          type="residual", sls=.05, aics=0,
                          pr=FALSE, what="observed-predicted",
                          tol=1e-12, maxiter=15, rel.tolerance=1e-5,
                          maxdim=5, ...)
{
  call <- match.call()
  cmethod <- match.arg(cmethod)
  if(cmethod=='hare')
    {
      require('polspline') ||
      {
        cat('polspline package not installed.  Reverting to cmethod="KM"\n')
        cmethod <- 'KM'
      }
    }
  
  if(!length(fit$y)) stop("fit did not store y")
  oldopt <- options(digits=3)
  on.exit(options(oldopt))
  unit <- fit$units
  if(unit=="") unit <- "Day"
  ny      <- dim(fit$y)
  nevents <- sum(fit$y[, ny[2]])

  ##Note: fit$y already has been transformed by the link function by psm

  survival <- survest.psm(fit, times=u, conf.int=FALSE)$surv

  if(cmethod=='hare' && missing(pred))
    {
      lim <- datadist(survival)$limits[c('Low:prediction','High:prediction'),]
      pred <- seq(lim[1], lim[2], length=100)
    }
  if(cmethod=='KM' && missing(cuts))
    {
      g <- max(1, floor(ny[1]/m))
      cuts <- quantile(c(0, 1, survival), seq(0, 1, length=g+1), na.rm=TRUE)
    }
  if(cmethod=='hare') cuts <- NULL
  else
    pred <- NULL

  dist <- fit$dist
  inverse <- survreg.distributions[[dist]]$itrans
  if(!length(inverse)) inverse <- function(x) x
  parms <- fit$parms
  
  distance <- function(x, y, fit, iter, u, fit.orig, what="observed", inverse,
                       pred, orig.cuts, maxdim, ...)
    {
      ##Assumes y is matrix with 1st col=time, 2nd=event indicator
      if(sum(y[,2]) < 5) return(NA)
      oldClass(fit) <- 'psm'   # for survest.psm which uses Survival.psm
      fit$dist <- fit.orig$dist
      psurv <- survest.psm(fit, linear.predictors=x,
                           times=u, conf.int=FALSE)$surv
      ##Assumes x really= x * beta
      if(length(orig.cuts))
        {
          pred.obs <- 
            groupkm(psurv, Surv(inverse(y[,1]), y[,2]), u=u, cuts=orig.cuts)
          dist <- if(what=="observed") pred.obs[,"KM"]
          else                         pred.obs[,"KM"] - pred.obs[,"x"]
        }
      else
        {
          pred.obs <- val.surv(fit, S=Surv(inverse(y[,1]), y[,2]), u=u,
                               est.surv=psurv,
                               pred=pred, maxdim=maxdim)
          dist <- if(what=='observed') pred.obs$actualseq
          else                         pred.obs$actualseq - pred
        }
      
      if(iter==0) storeTemp(pred.obs)
      dist
    }

  b <- min(10, B)
  overall.reps <- max(1, round(B/b))
  ## Bug in S prevents>10 loops in predab.resample
  if(pr)
    cat("\nAveraging ", overall.reps," repetitions of B=",b,"\n\n")
  rel  <- 0
  opt  <- 0
  nrel <- 0
  B    <- 0
  
  for(i in 1:overall.reps)
    {
      reliability <-
        predab.resample(fit, method=method,
                        fit=survreg.fit2, measure=distance,
                        pr=pr, B=b, bw=bw, rule=rule, type=type,  
                        u=u, m=m, what=what,
                        dist=dist, inverse=inverse, parms=parms,
                        fixed=fixed, family=family,
                        sls=sls, aics=aics, strata=FALSE,
                        tol=tol, pred=pred, orig.cuts=cuts, maxiter=maxiter,
                        rel.tolerance=rel.tolerance, maxdim=maxdim, ...)
      n <- reliability[,"n"]
      rel <- rel + n * reliability[,"index.corrected"]
      opt <- opt + n * reliability[,"optimism"]
      nrel <- nrel + n
      B <- B + max(n)	
      if(pr) print(reliability)
    }

  mean.corrected <- rel/nrel
  mean.opt       <- opt/nrel
  rel <- cbind(mean.optimism=mean.opt, mean.corrected=mean.corrected, n=nrel)
  if(pr)
    {
      cat("\nMean over ",overall.reps," overall replications\n\n")
      print(rel)
    }

  if(cmethod=='KM')
    {
      pred <- pred.obs[,"x"]
      KM   <- pred.obs[,"KM"]
      se   <- pred.obs[,"std.err"]
      obs.corrected <- KM - mean.opt
      structure(cbind(reliability[,c("index.orig","training","test"),
                                  drop=FALSE],
                      rel,mean.predicted=pred, KM=KM,
                      KM.corrected=obs.corrected, std.err=se),
                predicted=survival,
                class="calibrate", u=u, units=unit, n=ny[1], d=nevents, 
                p=length(fit$coefficients)-1, m=m, B=B, what=what, call=call)
    }
  else
    {
      calibrated            <- pred.obs$actualseq
      calibrated.corrected  <- calibrated - mean.opt
      structure(cbind(pred=pred,
                      reliability[,c("index.orig","training","test"),
                                  drop=FALSE],
                      rel, calibrated=calibrated,
                      calibrated.corrected=calibrated.corrected),
                predicted=survival,
                class="calibrate", u=u, units=unit, n=ny[1], d=nevents, 
                p=length(fit$coefficients)-1, m=m, B=B, what=what, call=call)
    }
}
