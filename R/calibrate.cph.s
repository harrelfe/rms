#Resampling optimism of reliability of a Cox survival model
#For predicting survival at a fixed time u, getting grouped K-M estimates
#with avg. of m subjects in a group, or using cutpoints cuts if present,
#e.g. cuts=c(0,.2,.4,.6,.8,1).
#B: # reps  method=see predab.resample
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each rep
#what="observed" to get optimism in observed (Kaplan-Meier) survival for
#groups by predicted survival
#what="observed-predicted" to get optimism in KM - Cox - more suitable if
#distributions of predicted survival vary greatly withing quantile groups
#defined from original sample's predicted survival

calibrate.cph <- function(fit, cmethod=c('hare', 'KM'),
                          method="boot", u, m=150, pred, cuts, B=40,
                          bw=FALSE, rule="aic",
                          type="residual", sls=.05, aics=0, force=NULL,
                          estimates=TRUE, pr=FALSE, what="observed-predicted",
                          tol=1e-12, maxdim=5, ...)
{
  call    <- match.call()
  cmethod <- match.arg(cmethod)
  ## if(cmethod=='hare')
  ##  {
  ##    require('polspline') ||
  ##    {
  ##      cat('polspline package not installed.  Reverting to cmethod="KM"\n')
  ##      cmethod <- 'KM'
  ##    }
  ##  }

  oldopt <- options('digits')
  options(digits=3)
  on.exit(options(oldopt))
  unit <- fit$units
  if(unit=="") unit <- "Day"
  
  ssum <- fit$surv.summary
  if(!length(ssum)) stop('did not use surv=TRUE for cph( )')
  cat("Using Cox survival estimates at ", dimnames(ssum)[[1]][2],
      " ", unit, "s\n", sep="")
  
  surv.by.strata <- ssum[2, , 1] #2nd time= at u, all strata
  xb <- fit$linear.predictors
  if(length(stra <- fit$strata)) 
    surv.by.strata <- surv.by.strata[stra]
  survival <- as.vector(surv.by.strata ^ exp(xb))

  if(cmethod=='hare' && missing(pred)) {
    lim <- datadist(survival)$limits[c('Low:prediction','High:prediction'),]
    pred <- seq(lim[1], lim[2], length=100)
  }
  if(cmethod=='KM' && missing(cuts)) {
    g <- max(1, floor(length(xb) / m))
    cuts <- unique(quantile(c(0, 1, survival), seq(0, 1, length=g + 1),
                            na.rm=TRUE))
  }
  if(cmethod=='hare') cuts <- NULL
  else pred <- NULL
  
  distance <- function(x, y, strata, fit, iter, u, fit.orig, what="observed",
                       pred, orig.cuts, maxdim, ...) {
    ## Assumes y is matrix with 1st col=time, 2nd=event indicator

    if(sum(y[, 2]) < 5) return(NA)
    surv.by.strata <- fit$surv.summary[2, , 1]
    ##2 means to use estimate at first time past t=0 (i.e., at u)
    
    if(length(strata))
      surv.by.strata <- surv.by.strata[strata] #Get for each stratum in data
    cox <- as.vector(surv.by.strata ^ exp(x - fit$center))
    ##Assumes x really= x * beta

    if(length(orig.cuts)) {
      pred.obs <- groupkm(cox, Surv(y[,1], y[,2]), u=u, cuts=orig.cuts)
      dist <- if(what == "observed") pred.obs[, "KM"]
              else                   pred.obs[, "KM"] - pred.obs[, "x"]
    } else {
      pred.obs <- val.surv(fit, S=Surv(y[, 1], y[, 2]), u=u,
                           est.surv=cox,
                           pred=pred, maxdim=maxdim)
      dist <- if(what=='observed') pred.obs$actualseq
              else                 pred.obs$actualseq - pred
    }
    
    if(iter == 0 && pr) print(pred.obs)

    if(iter == 0) structure(dist, keepinfo=list(pred.obs=pred.obs)) else
    dist
  }
  
  coxfit <- function(x, y, strata, u, iter=0, ...) {
    etime <- y[,1]
    e     <- y[,2]
    
    if(sum(e) < 5) return(list(fail=TRUE))
    x <- x	#Get around lazy evaluation creating complex expression
    f <- if(length(x)) {
      if(length(strata))
           cph(Surv(etime,e) ~ x + strat(strata), surv=TRUE, time.inc=u)
      else cph(Surv(etime,e) ~ x,                 surv=TRUE, time.inc=u)
    }
      else cph(Surv(etime,e) ~ strat(strata),     surv=TRUE, time.inc=u)
    ## Gets predicted survival at times 0, u, 2u, 3u, ...
    
    attr(f$terms, "Design") <- NULL
    ## Don't fool fastbw called from predab.resample
    f
  }

  b <- min(10, B)
  overall.reps <- max(1, round(B / b))
  ## Bug in S prevents>10 loops in predab.resample
  if(pr) cat("\nAveraging ", overall.reps, " repetitions of B=", b, "\n\n")
  rel  <- 0
  opt  <- 0
  nrel <- 0
  B    <- 0

  for(i in 1 : overall.reps) {
    reliability <-
      predab.resample(fit, method=method,
                      fit=coxfit, measure=distance,
                      pr=pr, B=b, bw=bw, rule=rule, type=type,  
                      u=u, m=m, what=what, sls=sls, aics=aics,
                      force=force, estimates=estimates,
                      pred=pred, orig.cuts=cuts, tol=tol, maxdim=maxdim, ...)
    kept     <- attr(reliability, 'kept') # TODO: accumulate over reps
    keepinfo <- attr(reliability, 'keepinfo')
    n    <- reliability[, "n"]
    rel  <- rel  + n * reliability[, "index.corrected"]
    opt  <- opt  + n * reliability[, "optimism"]
    nrel <- nrel + n
    B    <- B    + max(n)	
  }

  mean.corrected <- rel / nrel
  mean.opt       <- opt / nrel
  rel <- cbind(mean.optimism=mean.opt, mean.corrected=mean.corrected, n=nrel)
  if(pr) {
    cat("\nMean over ", overall.reps, " overall replications\n\n")
    print(rel)
  }
  
  e <- fit$y[, 2]
  pred.obs <- keepinfo$pred.obs
  if(cmethod == 'KM') {
    mean.predicted <- pred.obs[,"x"]
    KM             <- pred.obs[,"KM"]
    obs.corrected  <- KM - mean.opt
    
    structure(cbind(reliability[,c("index.orig","training","test"),
                                drop=FALSE],
                    rel, mean.predicted=mean.predicted, KM=KM,
                    KM.corrected=obs.corrected,
                    std.err=pred.obs[, "std.err", drop=FALSE]),
              predicted=survival, kept=kept,
              class="calibrate", u=u, units=unit, n=length(e), d=sum(e),
              p=length(fit$coefficients), m=m, B=B, what=what, call=call)
  } else {
    calibrated            <- pred.obs$actualseq
    calibrated.corrected  <- calibrated - mean.opt
    
    structure(cbind(pred=pred,
                    reliability[, c("index.orig", "training", "test"),
                                drop=FALSE],
                    rel, calibrated=calibrated,
                    calibrated.corrected=calibrated.corrected),
              predicted=survival, kept=kept,
              class="calibrate", u=u, units=unit, n=length(e), d=sum(e),
              p=length(fit$coefficients), m=m, B=B, what=what, call=call)
  }
}
