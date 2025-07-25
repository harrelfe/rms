predab.resample <-
  function(fit.orig,
           fit,
           measure, 
           method=c("boot", "crossvalidation", ".632", "randomization"),
           bw=FALSE,
           B=50,
           conf.int=0.95,
           pr=FALSE, prmodsel=TRUE,
           rule="aic",
           type="residual",
           sls=.05,
           aics=0,
           tol=.Machine$double.eps,
           force=NULL,
           estimates=TRUE,
           non.slopes.in.x=TRUE,
           kint=1,
           cluster,
           subset,
           group=NULL,
           allow.varying.intercepts=FALSE,
           debug=FALSE,
           saveraw=FALSE,
           ...)
{
  method <- match.arg(method)
  oldopt <- options('digits')
  options(digits=4)
  on.exit(options(oldopt))

  efit <- function(...) {
    message(...)
    list(fail=TRUE)
  }

  if(getOption('rmsdebug', FALSE)) tryCatch <- function(x, ...) x

  if(saveraw) genv <- .GlobalEnv
  
  ## Following logic prevents having to load a copy of a large x object
  if(any(match(c("x", "y"), names(fit.orig), 0) == 0))
    stop("must have specified x=T and y=T on original fit")
  
  fparms  <- fit.orig[c("terms", "Design")]
  oassign <- fit.orig$assign
  
  non.slopes <- num.intercepts(fit.orig, 'coef')

#  x.index <- if(non.slopes==0 || non.slopes.in.x) function(i,...) i
#  else
#    function(i, ns) {
#      if(any(i > ns)) i[i > ns] - ns
#      else NULL
#    }
  x.index <- function(i, ns)
    if(ns == 0) i else setdiff(i, 1 : ns) - ns
  
  Xb <- function(x, b, non.slopes, n, kint=1) {
    if(length(x)) matxv(x, b, kint=kint)
    else if(non.slopes == 0 || ! length(kint)) rep(0, n)
    else rep(b[kint], n)
  }
#    if(length(x)) {
#      if(non.slopes == 0 || non.slopes.in.x) x %*% b
#      else b[kint] + x %*% b[-(1 : non.slopes)]
#    }
#    else {
#      if(non.slopes==0) rep(0, n)
#      else
#        rep(b[kint], n)
#    }
#  }

  nac <- fit.orig$na.action
  
  x <- as.matrix(fit.orig[['x']])
  n <- nrow(x)
  
  ## Remove model.matrix class for subset operations later
  attr(x,'class') <- NULL	

  y <- fit.orig[['y']]
  if(! inherits(y, 'Ocens') && ! is.Surv(y)) {
    # if(! is.factor(y)) y <- factor(y)            ## ?? why was this ever here?
    # y <- as.matrix(unclass(y) - 1L)
    y <- as.matrix(unclass(y))
  }
  
  ## some subjects have multiple records now
  multi <- ! missing(cluster)

  if(length(group)) {
    if(multi || method != 'boot')
      stop('group is currently allowed only when method="boot" and cluster is not given')
    
    if(length(group) > n) {
      ## Missing observations were deleted during fit
      if(length(nac)) 
        j <- ! is.na(naresid(nac, y) %*% rep(1, ncol(y)))
      group <- group[j]
    }
    
    if(length(group) != n)
      stop('length of group does not match # rows used in fit')
    
    group.inds <- split(1 : n, group)  # see bootstrap()
    ngroup <- length(group.inds)
  }
  else ngroup <- 0
  
  if(multi) {
    if(method != 'boot')
      stop('cluster only implemented for method="boot"')
    
    if(length(cluster) > n) {
      ## Missing observations were deleted during fit
      if(length(nac)) {
        j <- ! is.na(naresid(nac, y) %*% rep(1, ncol(y)))
        cluster <- cluster[j]
      }
    }

    if(length(cluster) != n)
      stop('length of cluster does not match # rows used in fit')
    
    if(any(is.na(cluster)))
      stop('cluster has NAs')
    
    n.orig <- length(unique(cluster))
    cl.samp <- split(1 : n, cluster)
  }
  else n.orig <- n

  if(! missing(subset)) {
    if(length(subset) > n && length(nac)) {
      j <- ! is.na(naresid(nac, y) %*% rep(1, ncol(y)))
      subset <- subset[j]
    }
    
    if(length(subset) != n  && all(subset >= 0))
      stop('length of subset does not match # rows used in fit')
    
    if(any(is.na(subset))) stop('subset has NAs')
    
    if(! is.logical(subset)) {
      subset2         <- rep(FALSE, n)
      subset2[subset] <- TRUE
      subset          <- subset2
      subset2         <- NULL
    }
  }
  
  stra <- fit.orig$strata

  if(bw) {
    if(fit.orig$fail) stop('Original fit failed')

    if(prmodsel) cat("\n		Backwards Step-down - Original Model\n")
    fbw <- fastbw(fit.orig, rule=rule, type=type, sls=sls, aics=aics,
                  eps=tol, force=force)

    if(prmodsel) print(fbw, estimates=estimates)

    orig.col.kept <- fbw$parms.kept
    if(! length(orig.col.kept))
      stop("no variables kept in original model")
    ## Check that x.index works if allow.varying.intercepts
    
    xcol <- x.index(orig.col.kept, non.slopes)
    ## Refit subset of predictors on whole sample
    fit.orig <- fit(x[, xcol, drop=FALSE], y, strata=stra,
                    iter=0, tol=tol, xcol=xcol, ...)
    if(length(fit.orig$fail) && fit.orig$fail)
      stop('Refitting the stepdown model on the whole sample failed')
    
  }
  else orig.col.kept <- seq(along=fit.orig$coef)

  b <- fit.orig$coef
  xcol <- x.index(orig.col.kept, non.slopes)
  xb <- Xb(x[, xcol, drop=FALSE], b, non.slopes, n,
           kint=kint)

  index.orig <- if(missing(subset))
    measure(xb, y, strata=stra, fit=fit.orig, iter=0, evalfit=TRUE,
            fit.orig=fit.orig,
            kint=kint, ...)
  else 
    measure(xb[subset], y[subset,,drop=FALSE], strata=stra[subset],
            fit=fit.orig,
            iter=0, evalfit=FALSE, fit.orig=fit.orig, kint=kint, ...)
  keepinfo <- attr(index.orig, 'keepinfo')
  
  nindex     <- length(index.orig)
  test.stat  <- double(nindex)
  train.stat <- test.stat
  btrain     <- btest <- matrix(NA, B, nindex)
  name       <- fparms$Design$name
  if(bw) varin <- matrix(FALSE, nrow=B, ncol=length(name))

  # Assuming a naming convention, find out which elements of index.orig correspond to
  # calibration intercept or slope
  # No names on index.orig if validating nonparametric calibration curves

#  index_type <- if(! length(names(index.orig))) rep('npcal', length(index.orig))
#  else ifelse(tolower(names(index.orig)) %in% c('intercept', 'slope', 'emax'), 'lcal', 'other')

  # index_type  Meaning
  # npcal       points on nonparametric calibration curve
  # lcal        linear calibration parameters (slope, intercept, Emax)
  # other       regular model performance metrics (rank correlation, Brier, etc.)
  
  j   <- 0
  num <- 0

  if(method == "crossvalidation") {
    per.group <- n / B
    if(per.group < 2) {
      stop("B > n/2")
    }
    
    sb <- sample(n, replace=FALSE)
  }
  ##Cross-val keeps using same random set of indexes, without replacement
  
  ntest <- 0 #Used in getting weighted average for .632 estimator

  if(method == ".632") {
    ## Must do assignments ahead of time so can weight estimates
    ## according to representation in bootstrap samples
    S <- matrix(integer(1), nrow=n, ncol=B)
    W <- matrix(TRUE, nrow=n, ncol=B)
    for(i in 1 : B) {
      S[, i] <- s <- sample(n, replace=TRUE)
      W[s, i] <- FALSE  #now these obs are NOT omitted
    }
    
    nomit <- drop(W %*% rep(1,ncol(W)))  #no. boot samples omitting each obs
    if(min(nomit) == 0)
      stop("not every observation omitted at least once ",
           "in bootstrap samples.\nRe--run with larger B")
    
    W <- apply(W / nomit, 2, sum) / n
    if(pr) {
      cat("\n\nWeights for .632 method (ordinary bootstrap weights ",
          format(1 / B), ")\n", sep="")
      print(summary(W))
    }
  }
  pb <- setPb(B, type=if(method == 'crossvalidation') 'Cross' else 'Boot',
              onlytk=! pr,
              every=1*(B < 20) + 5*(B >= 20 & B < 50) +
                10*(B >= 50 & B < 100) + 20*(B >= 100 & B < 1000) +
                50*(B >= 1000))
  for(i in 1 : B) {
    pb(i)
    
    switch(method,
           crossvalidation = {
             is <- 1 + round((i - 1) * per.group)
             ie <- min(n, round(is + per.group - 1))
             test <- sb[is : ie]
             train <- -test
           }, #cross-val
           boot = {
             if(ngroup) {
               train <- integer(n.orig)

               for(si in 1 : ngroup) {
                 gi <- group.inds[[si]]
                 lgi <- length(gi)
                 train[gi] <- if(lgi == 1) gi
                 else {
                   ## sample behaves differently when first arg is
                   ## a single integer
                   sample(gi, lgi, replace=TRUE)
                 }
               }
             }
             else {
               train <- sample(n.orig, replace=TRUE)
               if(multi) train <- unlist(cl.samp[train])
             }
             test <- 1 : n
           },    #boot
           ".632" = {
             train <- S[, i]
             test <- -train
           },   #boot .632
           randomization =
           {
             train <- sample(n, replace=FALSE)
             test <- 1 : n
           }
           )   #randomization

    xtrain <- if(method == "randomization") 1 : n
    else train

    if(debug) {
      cat('\nSubscripts of training sample:\n')
      print(train)
      cat('\nSubscripts of test sample:\n')
      print(test)
    }
    f <- tryCatch(fit(x[xtrain, , drop=FALSE], y[train, , drop=FALSE],
                      strata=stra[train], iter=i, tol=tol, ...), error=efit)
    if(! length(f$fail)) f$fail <- FALSE
    f$assign <- NULL  #Some programs put a NULL assign (e.g. ols.val fit)
    ni <- num.intercepts(f)
    
    fail <- f$fail
    if(! fail) {
      ## Following if..stop was before f$assign above
      if(! allow.varying.intercepts && ni != non.slopes) {
        stop('A training sample has a different number of intercepts (', ni ,')\n',
             'than the original model fit (', non.slopes, ').\n',
             'You probably fit an ordinal model with sparse cells and a re-sample\n',
             'did not select at least one observation for each value of Y.\n',
             'Add the argument group=y where y is the response variable.\n',
             'This will force balanced sampling on levels of y.')
      }
      
      clf <- attr(f, "class")  # class is removed by c() below

      f[names(fparms)] <- fparms
      assign <- oassign
      ## Slopes are shifted to the left when fewer unique values of Y
      ## occur (especially for orm models) resulting in fewer intercepts
      if(non.slopes != ni) for(z in 1 : length(assign))
        assign[[z]] <- assign[[z]] - (non.slopes - ni)
      f$assign <- assign
      attr(f, "class") <- clf
      if(! bw) {
        coef <- f$coef
        col.kept <- seq(along=coef)
      }
      else {
        f <- fastbw(f, rule=rule, type=type, sls=sls, aics=aics,
                    eps=tol, force=force)
        
        if(pr && prmodsel) print(f, estimates=estimates)
        
        varin[j + 1, f$factors.kept] <- TRUE
        col.kept <- f$parms.kept

        if(! length(col.kept))
          f <- tryCatch(fit(NULL, y[train,, drop=FALSE], stra=stra[xtrain],
                            iter=i, tol=tol,...), error=efit)
        else {
          xcol <- x.index(col.kept, ni)
          f <- tryCatch(fit(x[xtrain, xcol, drop=FALSE], strata=stra[xtrain],
                            y[train,, drop=FALSE],
                            iter=i, tol=tol, xcol=xcol, ...), error=efit)
        }

        if(length(f$fail) && f$fail) fail <- TRUE
        else coef <- f$coef
      }
    }
    
    if(! fail) {
      j    <- j + 1
      xcol <- x.index(col.kept, ni)
      xb   <- Xb(x[,xcol,drop=FALSE], coef, ni, n, kint=kint)
      
      if(missing(subset)) {
        train.statj <-
          measure(xb[xtrain], y[train,,drop=FALSE], strata=stra[xtrain], 
                  fit=f, iter=i, fit.orig=fit.orig, evalfit=TRUE,
                  kint=kint, ...)
        test.statj <- measure(xb[test], y[test,,drop=FALSE],
                              strata=stra[test],
                              fit=f, iter=i, fit.orig=fit.orig,
                              evalfit=FALSE,
                              kint=kint, ...)
      }
      else {
        ii <- xtrain
        
        if(any(ii < 0)) ii <- (1 : n)[ii]
        
        ii <- ii[subset[ii]]
        train.statj <- measure(xb[ii], y[ii,,drop=FALSE],
                               strata=stra[ii],
                               fit=f, iter=i, fit.orig=fit.orig,
                               evalfit=FALSE,
                               kint=kint, ...)
        
        ii <- test
        if(any(ii < 0)) ii <- (1 : n)[ii]
        
        ii <- ii[subset[ii]]
        test.statj <- measure(xb[ii], y[ii,,drop=FALSE], fit=f,
                              iter=i, strata=stra[ii],
                              fit.orig=fit.orig, evalfit=FALSE,
                              kint=kint, ...)
      }
      
#      stt <- ifelse(index_type %in% c('lcal', 'npcal'), test.statj, train.statj)
      btrain[j, ] <- train.statj
      btest [j, ] <- test.statj

      # Type of Index  Bootstrap Estimates Computed for CLs
      # lcal           estimates from test (original) sample
      # npcal          " " "
      # other          estimates from training (bootstrap) sample

      na <- is.na(train.statj + test.statj)
      num <- num + ! na
      if(pr) 
        print(cbind(training=train.statj, test=test.statj))
      
      # For a given resample, neither training nor test estimates are counted
      # if either one is NA
      train.statj[na] <- 0
      test.statj[na]  <- 0
      if(method == ".632") {
        wt <- W[i]
        if(any(na))
          warning('method=".632" does not properly handle missing summary indexes')
      }
      else wt <- 1
      
      train.stat <- train.stat + train.statj
      test.stat  <- test.stat  + test.statj * wt
      ntest      <- ntest + 1
    } 
  }
  
  if(pr) cat("\n\n")
  
  if(pr && (j != B))
    cat("\nDivergence or singularity in", B - j, "samples\n")

  btrain <- btrain[1 : j,, drop=FALSE]
  btest  <- btest [1 : j,, drop=FALSE]
  
  train.stat <- train.stat / num
  
  if(method != ".632") {
    test.stat <- test.stat / num
    optimism <- train.stat - test.stat
  }
  else optimism <- .632 * (index.orig - test.stat)
  
  res <- cbind(index.orig=index.orig, training=train.stat, test=test.stat,
               optimism=optimism, index.corrected=index.orig - optimism)

  if(saveraw)
    assign('.predab_raw.', list(orig=index.orig, btrain=btrain, btest=btest), envir=genv)
  
  if(conf.int > 0) {
    z  <- qnorm((conf.int + 1) / 2)
    s <- apply(btrain - 1.25 * btest, 2, dualSD)   # dualSD is in Hmisc
    # Reverse top and bottom since the optimism is subtracted from the apparent index
    Lower <- index.orig - optimism - z * s['top',   ]
    Upper <- index.orig - optimism + z * s['bottom',]
    res   <- cbind(res, Lower, Upper)
  }
  res <- cbind(res, n=num)

  # Note: For npcal, Lower & Upper are on calibrated probabilities minus
  # predicted probabilities (predy); predy needs to be added back in
  # calibrate* functions
  
  if(bw) {
    varin <- varin[1 : j, , drop=FALSE]
    dimnames(varin) <- list(rep("", j), name)
  }
  structure(res, class='validate', kept=if(bw) varin, keepinfo=keepinfo)
}

