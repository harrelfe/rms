## This is a modification of the R survival package's coxph function
## written by Terry Therneau and ported to R by Thomas Lumley
cph <- function(formula     = formula(data),
                data        = environment(formula),
                weights,
                subset,
                na.action   = na.delete, 
                method      =c("efron", "breslow", "exact",
                               "model.frame", "model.matrix"),
                singular.ok = FALSE,
                robust      = FALSE,
                model       = FALSE,
                x           = FALSE,
                y           = FALSE,
                se.fit      = FALSE,
                linear.predictors = TRUE,
                residuals         = TRUE,
                nonames           = FALSE,
                eps         = 1e-4,
                init,
                iter.max    = 10,
                tol         = 1e-9,
                surv        = FALSE,
                time.inc,
                type        = NULL,
                vartype     = NULL,
                debug       = FALSE,
                ...)
{
  method <- match.arg(method)
  call   <- match.call()

  if (! inherits(formula,"formula")) {
    ## I allow a formula with no right hand side
    ## The dummy function stops an annoying warning message "Looking for
    ##  'formula' of mode function, ignored one of mode ..."
    if (inherits(formula, "Surv")) {
      xx <- function(x) formula(x)
      formula <- xx(paste(deparse(substitute(formula)), 1, sep="~"))
    }
    else stop("Invalid formula")
  }

  data <-
    modelData(data, formula,
              weights=if(! missing(weights)) eval(substitute(weights), data),
              subset =if(! missing(subset )) eval(substitute(subset),  data),
              na.action=na.action, dotexpand=FALSE)

  nstrata <- 0
  Strata <- NULL

  odb <- .Options$debug
  if(length(odb) && is.logical(odb) && odb) debug <- TRUE

  if(length(z <- attr(terms(formula, allowDotAsName=TRUE), "term.labels")) > 0
     && any(z !=".")) { #X's present

    X    <- Design(data, formula, specials=c('strat', 'strata'))
    
    atrx       <- attributes(X)
    atr        <- atrx$Design
    nact       <- atrx$na.action
    sformula   <- atrx$sformula
    mmcolnames <- atr$mmcolnames
    if(method == "model.frame") return(X)

    Terms <- terms(sformula, specials=c('strat', 'strata'), data=data)
  

    asm   <- atr$assume.code
    name  <- atr$name

    specials <- attr(Terms, 'specials')
    if(length(specials$strata)) stop('cph supports strat(), not strata()')
    stra    <- specials$strat

    cluster <- attr(X, 'cluster')
    if(length(cluster)) {
      if(missing(robust)) robust <- TRUE
      attr(X, 'cluster') <- NULL
    }
    
    Terms.ns     <- Terms
    if(length(stra)) {
      temp <- untangle.specials(Terms.ns, "strat", 1)
      Terms.ns <- Terms.ns[- temp$terms]	#uses [.terms function
      
      Strata <- list()
      strataname <- attr(Terms, 'term.labels')[stra - 1]
      
      j <- 0
      for(i in (1 : length(asm))[asm == 8]) {
        nstrata <- nstrata + 1
        xi <- X[[i + 1]]
        levels(xi) <- paste(name[i], "=", levels(xi), sep="")
        Strata[[nstrata]] <- xi
      }
      Strata <- interaction(as.data.frame(Strata), drop=TRUE)
    }
    xpres <- length(asm) && any(asm != 8)
    Y <- model.extract(X, 'response')
    if(! inherits(Y, "Surv"))
      stop("response variable should be a Surv object")
    n <- nrow(Y)
    weights <- model.extract(X, 'weights')
    offset  <- attr(X, 'offset')
    ##  Cox ph fitter routines expect null if no offset
    
    ##No mf if only strata factors
    if(! xpres) {
      X <- matrix(nrow=0, ncol=0)
      assign <- NULL
    }
    else {
      X <- model.matrix(sformula, X)
      ## Handle special case where model was fitted using previous fit$x
      alt <- attr(mmcolnames, 'alt')
      if(debug) {
        prn(sformula)
        print(cbind('colnames(X)'=colnames(X)[-1],
                    mmcolnames=mmcolnames,
                    'Design colnames'=atr$colnames,
                    alt=alt))
      }
#        prn(colnames(X)); prn(mmcolnames); prn(alt)}
      if(! all(mmcolnames %in% colnames(X)) && length(alt)) mmcolnames <- alt
      X <- X[, mmcolnames, drop=FALSE]
      assign <- attr(X, "assign")
      assign[[1]] <- NULL  # remove intercept position, renumber
    }
    
    nullmod <- FALSE
  }
  else {	## model with no right-hand side
    X        <- NULL
    Y        <- data[[1]]
    sformula <- formula
    mmcolnames <- ''
    weights  <- if('(weights)' %in% names(data)) data[['(weights)']]
    atr <- atrx  <- NULL
    Terms    <- terms(formula, allowDotAsName=TRUE)
    if(! inherits(Y, "Surv"))
      stop("response variable should be a Surv object")
    
    Y <- Y[! is.na(Y)]
    assign  <- NULL
    xpres   <- FALSE
    nullmod <- TRUE
    nact    <- NULL
  }
  ny <- ncol(Y)
  maxtime <- max(Y[, ny - 1])

  rnam <- if(! nonames) dimnames(Y)[[1]]
  if(xpres) dimnames(X) <- list(rnam, atr$colnames)

  if(method == "model.matrix") return(X)

  time.units <- units(Y)
  if(! length(time.units) || time.units == '') time.units <- "Day"
  
  if(missing(time.inc)) {
    time.inc <- switch(time.units,
                       Day   = 30,
                       Month = 1,
                       Year  = 1,
                       maxtime / 10)
    
    if(time.inc >= maxtime | maxtime / time.inc > 25)
      time.inc <- max(pretty(c(0, maxtime))) / 10
  }

  ytype <- attr(Y, 'type')
  if(nullmod) f <- NULL
  else {
    fitter <-
      if( method == "breslow" || method  == "efron") {
        if (ytype ==  'right') coxph.fit
        else agreg.fit
      }
      else if (method == 'exact') {
        if(ytype == 'right') getFromNamespace('coxexact.fit', 'survival')
        else
          agexact.fit
        }
      else
        stop(paste ("Unknown method", method))
    
    if (missing(init)) init <- NULL

    f <- fitter(X, Y,
                strata=Strata, offset=offset,
                weights=weights, init=init,
                method=method, rownames=rnam,
                control=coxph.control(eps=eps, toler.chol=tol,
                                      toler.inf=1, iter.max=iter.max))
  }
  if (is.character(f)) {
    cat("Failure in cph:\n", f, "\n")
    return(structure(list(fail=TRUE), class="cph"))
  }
  else {
    if(length(f$coefficients) && any(is.na(f$coefficients))) {
      vars <- names(f$coefficients)[is.na(f$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable",
                   paste(vars, collapse=" "))
      if(singular.ok) warning(msg)
      else {
        cat(msg,"\n")
        return(structure(list(fail=TRUE), class="cph"))
      }
    }
  }
  f$terms      <- Terms
  f$sformula   <- sformula
  f$mmcolnames <- mmcolnames
  
  if(robust) {
    f$naive.var <- f$var
    ## Terry gets a little tricky here, calling resid before adding
    ## na.action method to avoid re-inserting NAs.  Also makes sure
    ## X and Y are there
    if(! length(cluster)) cluster <- FALSE
    
    fit2 <- c(f, list(x=X, y=Y, weights=weights, method=method))
    if(length(stra)) fit2$strata <- Strata
    
    r <- getS3method('residuals', 'coxph')(fit2, type='dfbeta',
                                           collapse=cluster, weighted=TRUE)
    f$var <- t(r) %*% r
  }
  
  nvar <- length(f$coefficients)

  ev <- factor(Y[, ny], levels=0 : 1, labels=c("No Event", "Event"))
  n.table <- {
    if(! length(Strata)) table(ev, dnn='Status')
    else table(Strata, ev, dnn=c('Stratum', 'Status'))
  }
  f$n <- n.table
  nevent <- sum(Y[, ny])
  if(xpres) {
    logtest <- -2 * (f$loglik[1] - f$loglik[2])
    R2.max  <-  1 - exp(2 * f$loglik[1] / n)
    R2 <- (1 - exp(- logtest / n)) / R2.max
    P  <- 1 - pchisq(logtest,nvar)
    gindex <- GiniMd(f$linear.predictors)
    dxy <- dxy.cens(f$linear.predictors, Y, type='hazard')['Dxy']
    stats <- c(n, nevent, logtest, nvar, P, f$score, 
               1-pchisq(f$score,nvar), R2, dxy, gindex, exp(gindex))
    names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P", 
                      "Score", "Score P", "R2", "Dxy", "g", "gr")
  }
  else {
    stats <- c(n, nevent)
    names(stats) <- c("Obs", "Events")
  }
  
  f$method <- NULL
  if(xpres)
    dimnames(f$var) <- list(atr$colnames, atr$colnames)
  
  f <- c(f, list(call=call, Design=atr,
                 assign=DesignAssign(atr, 0, atrx$terms),
                 na.action=nact,
                 fail = FALSE, non.slopes = 0, stats = stats, method=method,
                 maxtime = maxtime, time.inc = time.inc,
                 units = time.units))

  if(xpres) {
    f$center <- sum(f$means * f$coefficients)
    f$scale.pred <- c("log Relative Hazard", "Hazard Ratio")
    attr(f$linear.predictors,"strata") <- Strata
    names(f$linear.predictors) <- rnam
    if(se.fit) {
      XX <- X - rep(f$means, rep.int(n, nvar))   # see scale() function
      ##  XX <- sweep(X, 2, f$means)	# center   (slower;so is scale)
      se.fit <- drop(((XX %*% f$var) * XX) %*% rep(1,ncol(XX)))^.5
      names(se.fit) <- rnam
      f$se.fit <- se.fit  	
    }
  }
  if(model) f$model <- data
  
  if(is.character(surv) || surv) {
    if(length(Strata)) {
      iStrata <- as.character(Strata)
      slev <- levels(Strata)
      nstr <- length(slev)
    } else nstr <- 1
    srv  <- NULL
    tim  <- NULL
    s.e. <- NULL
    timepts <- seq(0, maxtime, by=time.inc)
    s.sum <- array(double(1),
                   c(length(timepts), nstr, 3),
                   list(format(timepts), paste("Stratum", 1 : nstr),
                        c("Survival", "n.risk", "std.err")))
    g <- list(n=sum(f$n),
              coefficients=f$coefficients,
              linear.predictors=f$linear.predictors,
              method=f$method, type=type, means=f$means, var=f$var,
              x=X, y=Y, strata=Strata, offset=offset, weights=weights,
              terms=Terms, call=call)
    g <- survfit.cph(g, se.fit=is.character(surv) || surv,
                     type=type, vartype=vartype, conf.type='log')

    strt <- if(nstr > 1) rep(names(g$strata), g$strata)

    for(k in 1 : nstr) {
      j    <- if(nstr == 1) TRUE else strt == slev[k]
      yy   <- Y[if(nstr == 1) TRUE else iStrata == slev[k], ny - 1]
      maxt <- max(yy)
      ##n.risk from surv.fit does not have usual meaning if not Kaplan-Meier
      
      tt <- c(0,  g$time[j])
      su <- c(1,  g$surv[j])
      se <- c(NA, g$std.err[j])
      
      if(maxt > tt[length(tt)]) {
        tt <- c(tt, maxt)
        su <- c(su, su[length(su)])
        se <- c(se, NA)
      }

      kk <- 0
      for(tp in timepts) {
        kk <- kk + 1
        t.choice <- max((1 : length(tt))[tt <= tp+1e-6])
        if(tp > max(tt) + 1e-6 & su[length(su)] > 0) {
          Su <- NA
          Se <- NA
        }
        else {
          Su <- su[t.choice]
          Se <- se[t.choice]
        }
        
        n.risk <- sum(yy >= tp)
        s.sum[kk, k, 1 : 3] <- c(Su, n.risk, Se)
      }
      
      if(! is.character(surv)) {
        if(nstr == 1) {
          tim  <- tt
          srv  <- su
          s.e. <- se
        }
        else {
          tim  <- c(tim,  list(tt))
          srv  <- c(srv,  list(su))
          s.e. <- c(s.e., list(se))
        }
      }
    }
    
    if(is.character(surv)) f$surv.summary <- s.sum
    else {
      if(nstr > 1) {
        names(srv) <- names(tim) <- names(s.e.) <- levels(Strata) ###
      }
      
      f <- c(f, list(time=tim, surv=srv,
                     std.err=s.e., surv.summary=s.sum))		
    }
  }
  
  f$strata <- Strata    ### was $Strata
  if(x) f$x <- X
  if(y) f$y <- Y
  f$weights <- weights
  f$offset  <- offset

  if(! linear.predictors) f$linear.predictors <- NULL
  if(! residuals        ) f$residuals <- NULL
  class(f) <- c("cph", "rms", "coxph")
  f
}

coxphFit <- function(..., method, strata=NULL, rownames=NULL, offset=NULL,
                     init=NULL, toler.chol=1e-9, eps=.0001, iter.max=10,
                     type) {

  fitter <- if( method == "breslow" || method == "efron") {
              if (type == 'right') coxph.fit else agreg.fit
            }
            else if (method == 'exact') {
              if(type == 'right') getFromNamespace('coxexact.fit', 'survival')
              else
                agexact.fit
              }
  else stop("Unkown method ", method)

  res <- fitter(..., strata=strata, rownames=rownames,
                offset=offset, init=init, method=method,
                control=coxph.control(toler.chol=toler.chol, toler.inf=1,
                                      eps=eps, iter.max=iter.max))
  
  if(is.character(res)) return(list(fail=TRUE))
  
  if(iter.max > 1 && res$iter >= iter.max) return(list(fail=TRUE))

  res$fail <- FALSE
  res
}

Survival.cph <- function(object, ...) {
  if(! length(object$time) || ! length(object$surv))
    stop("did not specify surv=T with cph")
  f <- function(times, lp=0, stratum=1, type=c("step","polygon"),
                time, surv) {
    type <- match.arg(type)
    if(length(stratum) > 1) stop("does not handle vector stratum")
    if(length(times) == 0) {
      if(length(lp) > 1) stop("lp must be of length 1 if times=NULL")
      return(surv[[stratum]] ^ exp(lp))
    }
    s <- matrix(NA, nrow=length(lp), ncol=length(times),
                dimnames=list(names(lp), format(times)))
    if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
      if(type == "polygon") {
        if(length(lp) > 1 && length(times) > 1)
          stop('may not have length(lp)>1 & length(times>1) when type="polygon"')
        su <- approx(time, surv, times, ties=mean)$y
        return(su ^ exp(lp))
      }
    for(i in 1 : length(times)) {
      tm <- max((1 : length(time))[time <= times[i] + 1e-6])
      su <- surv[tm]
      if(times[i] > max(time) + 1e-6) su <- NA
      s[,i] <- su ^ exp(lp)
    }
    drop(s)
  }
  formals(f) <- list(times=NULL, lp=0, stratum=1,
                     type=c("step","polygon"),
                     time=object$time, surv=object$surv)
  f
}

Quantile.cph <- function(object, ...) {
  if(! length(object$time) || ! length(object$surv))
    stop("did not specify surv=T with cph")
  f <- function(q=.5, lp=0, stratum=1, type=c("step","polygon"), time, surv) {
    type <- match.arg(type)
    if(length(stratum)>1) stop("does not handle vector stratum")
    if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
    Q <- matrix(NA, nrow=length(lp), ncol=length(q),
                dimnames=list(names(lp), format(q)))
    for(j in 1 : length(lp)) {
      s <- surv^exp(lp[j])
      if(type == "polygon") Q[j,] <- approx(s, time, q, ties=mean)$y
      else for(i in 1 : length(q))
        if(any(s <= q[i])) Q[j,i] <- min(time[s <= q[i]])  #is NA if none
    }
    drop(Q)
  }
  formals(f) <- list(q=.5, lp=0, stratum=1,
                     type=c('step','polygon'),
                     time=object$time, surv=object$surv)
  f
}


Mean.cph <- function(object, method=c("exact","approximate"),
                     type=c("step","polygon"), n=75, tmax=NULL, ...) {
  method <- match.arg(method)
  type   <- match.arg(type)
  
  if(! length(object$time) || ! length(object$surv))
    stop("did not specify surv=TRUE with cph")
  
  if(method == "exact") {
    f <- function(lp=0, stratum=1, type=c("step","polygon"),
                  tmax=NULL, time, surv) {
      type <- match.arg(type)
      if(length(stratum) > 1) stop("does not handle vector stratum")
      if(is.list(time)) {time <- time[[stratum]]; surv <- surv[[stratum]]}
      Q <- lp
      if(! length(tmax)) {
        if(min(surv) > 1e-3)
          warning(paste("Computing mean when survival curve only defined down to",
                        format(min(surv)), "\n Mean is only a lower limit"))
        k <- rep(TRUE, length(time))
      }
      else {
        if(tmax > max(time)) stop(paste("tmax=", format(tmax),
                                      "> max follow-up time=",
                                      format(max(time))))
        k <- (1 : length(time))[time <= tmax]
      }
      for(j in 1 : length(lp)) {
        s <- surv ^ exp(lp[j])
        Q[j] <- if(type == "step") sum(c(diff(time[k]), 0) * s[k]) else 
          trap.rule(time[k], s[k])
      }
      Q
    }
    formals(f) <- list(lp=0, stratum=1,
                       type=c("step","polygon"),
                       tmax=tmax,
                       time=object$time, surv=object$surv)
    return(f)
  }
  else {
    lp     <- object$linear.predictors
    lp.seq <- if(length(lp)) lp.seq <- seq(min(lp), max(lp), length=n) else 0
    
    time   <- object$time
    surv   <- object$surv
    nstrat <- if(is.list(time)) length(time) else 1
    areas  <- list()
    
    for(is in 1 : nstrat) {
      tim <- if(nstrat == 1) time else time[[is]]
      srv <- if(nstrat == 1) surv else surv[[is]]
      if(! length(tmax)) {
        if(min(srv) > 1e-3)
          warning(paste("Computing mean when survival curve only defined down to",
                        format(min(srv)),
                        "\n Mean is only a lower limit"))
        k <- rep(TRUE, length(tim))
      }
          else {
            if(tmax > max(tim)) stop(paste("tmax=",format(tmax),
                                           "> max follow-up time=",
                                           format(max(tim))))
            k <- (1 : length(tim))[tim <= tmax]
          }
      ymean <- lp.seq
      for(j in 1 : length(lp.seq)) {
        s <- srv ^ exp(lp.seq[j])
        ymean[j] <- if(type == "step") sum(c(diff(tim[k]),0) * s[k]) else 
        trap.rule(tim[k], s[k])
      }
      areas[[is]] <- ymean
    }
    if(nstrat > 1) names(areas) <- names(time)

    ff <- function(lp=0, stratum=1, lp.seq, areas) {
      
      if(length(stratum) > 1) stop("does not handle vector stratum")
      area <- areas[[stratum]]
      if(length(lp.seq) == 1 && all(lp == lp.seq))
        ymean <- rep(area, length(lp))
      else ymean <- approx(lp.seq, area, xout=lp, ties=mean)$y
      if(any(is.na(ymean)))
        warning("means requested for linear predictor values outside range of linear\npredictor values in original fit")
      names(ymean) <- names(lp)
      ymean
    }
    formals(ff) <- list(lp=0, stratum=1, lp.seq=lp.seq, areas=areas)
  }
  ff
}

predict.cph <- function(object, newdata=NULL,
                        type=c("lp", "x", "data.frame", "terms", "cterms",
                          "ccterms", "adjto", "adjto.data.frame", "model.frame"),
                        se.fit=FALSE, conf.int=FALSE,
                        conf.type=c('mean','individual','simultaneous'),
                        kint=1,
                        na.action=na.keep, expand.na=TRUE,
                        center.terms=type=="terms", ...) {
  type <- match.arg(type)
  predictrms(object, newdata, type, se.fit, conf.int, conf.type,
             kint,
             na.action, expand.na, center.terms, ...)
}

print.cph <- function(x, digits=4, table=TRUE, conf.int=FALSE,
                      coefs=TRUE,
                      title='Cox Proportional Hazards Model', ...)
{ 
  k <- 0
  z <- list()
  
  if(length(zz <- x$na.action)) {
    k <- k + 1
    z[[k]] <- list(type=paste('naprint', class(zz)[1], sep='.'), list(zz))
  }
  
  if(table && length(x$n) && is.matrix(x$n)) {
    k <- k + 1
    z[[k]] <- list(type='print', list(x$n))
  }
  
  if(length(x$coef)) {
    stats <- x$stats
    ci <- x$clusterInfo
    misc <- reListclean(Obs   =stats['Obs'],
                     Events=stats['Events'],
                     'Cluster on' = ci$name,
                     Clusters = ci$n,
                     Center   = round(x$center, digits))
    lr   <- reListclean('LR chi2'     = stats['Model L.R.'],
                     'd.f.'        = stats['d.f.'],
                     'Pr(> chi2)'  = stats['P'],
                     'Score chi2'  = stats['Score'],
                     'Pr(> chi2)'  = stats['Score P'])
    disc <- reListclean(R2 = stats['R2'],
                     Dxy = stats['Dxy'],
                     g  = stats['g'],
                     gr = stats['gr'])
    k <- k + 1
    headings <- c('', 'Model Tests', 'Discrimination\nIndexes')
    data     <- list(misc,
                     c(lr,   c(2,NA,4,2,4)),
                     c(disc, 3))
    z[[k]] <- list(type='stats', list(headings=headings, data=data))
    
    beta <- x$coef
    se <- sqrt(diag(x$var))
    k <- k + 1
    z[[k]] <- list(type='coefmatrix',
                   list(coef = x$coef,
                        se   = sqrt(diag(x$var))))
    if(conf.int) {
          
      zcrit <- qnorm((1 + conf.int)/2)
      tmp <- cbind(exp(beta), exp( - beta), exp(beta - zcrit * se),
                   exp(beta + zcrit * se))
      dimnames(tmp) <- list(names(beta),
                            c("exp(coef)", "exp(-coef)",
                              paste("lower .",
                                    round(100 * conf.int, 2), sep = ""),
                              paste("upper .",
                                    round(100 * conf.int, 2), sep = "")))
      k <- k + 1
      z[[k]] <- list(type='print', list(tmp, digits=digits))
    }
  }
  
  prModFit(x, title=title,
           z, digits=digits, coefs=coefs, ...)
}
