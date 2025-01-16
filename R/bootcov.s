bootcov <- function(fit, cluster, B=200, fitter, coef.reps=TRUE,
                    loglik=FALSE, pr=FALSE,
                    group=NULL, stat=NULL, seed=sample(10000, 1), ytarget=NULL, ...) {

  coxcph <- inherits(fit,'coxph') || inherits(fit,'cph')

  nfit <- class(fit)[1]

  if(any(c('x', 'y') %nin% names(fit)))
    stop("you did not specify x=TRUE and y=TRUE in the fit")
  X <- fit$x
  Y <- fit$y

  if(length(stat) > 1) stop('stat may only contain one statistic name')

  sc.pres <- 'scale' %in% names(fit)
  ns      <- num.intercepts(fit)

  ## See if ordinal regression being done
  yu      <- fit$yunique
  ychar   <- is.character(yu)
  ordinal <- nfit == 'orm' || (nfit == 'lrm' && length(yu) > 2)

  if(length(ytarget) && nfit != 'orm') stop('ytarget applies only to orm fits'
  )
  if(nfit == 'orm' && length(ytarget)) {
    if(is.na(ytarget)) {
      iref    <- fit$interceptRef
      ytarget <- if(ychar) yu[-1][iref] else median(Y)
    } else {
      iref <- if(ychar) which(yu[-1] == ytarget) else which.min(abs(yu[-1] - ytarget))
      if(! length(iref)) stop('no intercept corresponds to ytarget=', ytarget)
    }
  }

  ## Someday need to add resampling of offsets, weights    TODO
  if(missing(fitter))
    fitter <- quickRefit(fit, what='fitter', ytarget=ytarget, storevals=FALSE)

  if(! length(fitter)) stop("fitter not valid")

  if(loglik) {
    oosl <- switch(nfit,
                   ols=oos.loglik.ols,
                   lrm=oos.loglik.lrm,
                   cph=oos.loglik.cph,
                   psm=oos.loglik.psm,
                   Glm=oos.loglik.Glm)

      if(!length(oosl))
        stop('loglik=TRUE but no oos.loglik method for model in rmsMisc')

      Loglik        <- double(B + 1)
      Loglik[B + 1] <- oosl(fit)
    }
  else Loglik <- NULL

  n        <- nrow(X)
  Cof      <- fit$coefficients
  intnames <- names(Cof)[1 : ns]
 
  if(nfit == 'orm' && length(ytarget)) {
    message('Keeping only intercept ', iref, ' (position for original sample) for ytarget=', ytarget)
    ikeep         <- c(iref, (ns + 1) : length(Cof))
    Cof           <- Cof[ikeep]
    names(Cof)[1] <- 'Intercept'
  }
  
  p     <- length(Cof)
  vname <- names(Cof)
  if(sc.pres) {
    p     <- p + 1L
    vname <- c(vname, "log scale")
  }

  bar <- rep(0, p)
  cov <- matrix(0, nrow=p, ncol=p, dimnames=list(vname,vname))
  if(coef.reps) coefs <- matrix(NA, nrow=B, ncol=p, dimnames=list(NULL, vname))
  if(length(stat)) stats <- numeric(B)
  nry <- rep(0, B)

  Y  <- as.matrix(if(is.factor(Y)) unclass(Y) else Y)
  ny <- ncol(Y)

  Strata <- fit$strata

  nac <- fit$na.action

  if(length(group)) {
    if(length(group) > n) {
      ## Missing observations were deleted during fit
      if(length(nac)) {
        j <- !is.na(naresid(nac, Y) %*% rep(1, ny))
        group <- group[j]
      }
    }

    if(length(group) != n)
      stop('length of group does not match # rows used in fit')

    group.inds <- split(1:n, group)
    ngroup     <- length(group.inds)
  }
  else ngroup <- 0

  # Given a vector of intercepts, with those corresponding to non-sampled y
  # equal to NA, use linear interpolation/extrapolation to fill in the NAs
  fillin <- function(y, alpha) {
    if(length(y) != length(alpha)) stop('lengths of y and alpha must match')
    i <- ! is.na(alpha)
    if(sum(i) < 2)
      stop('need at least 3 distinct Y values sampled to be able to extrapolate intercepts')
    est_alpha <- approxExtrap(y[i], alpha[i], xout=y[! i])$y
    alpha[! i] <- est_alpha
    alpha
  }

  process_ints <- function(cof, ints_fitted, nry) {
    if(nry == 0) return(cof)
    if(length(ytarget))
      stop('Program logic error: ytarget is specified but there are still ',
           'missing intercepts in a bootstrap sample')
    if((nfit == 'orm') && ! ychar) {
      # Numeric Y; use linear interpolation/extrapolation to fill in
      # intercepts for non-sampled Y values
      alphas <- structure(rep(NA, ns), names=intnames)  # intercept template
      ints_actual <- cof[1 : ints_fitted]
      alphas[names(ints_actual)] <- ints_actual
      if(sum(is.na(alphas)) != nry) stop('program logic error in alphas') 
      alphas <- fillin(yu[- 1], alphas)
      return(c(alphas, cof[- (1 : ints_fitted)]))
      }
    stop('Bootstrap sample did not include the following intercepts. ',
         'Do minimal grouping on y using ordGroupBoot() to ensure that ',
         'all bootstrap samples will have all original distinct y values ',
         'represented.  ', paste(setdiff(vname, names(cof)), collapse=' '))
    }

  set.seed(seed)

  if(missing(cluster)) {
    clusterInfo <- NULL
    nc <- n
    b  <- 0
    pb <- setPb(B, type='Boot', onlytk=! pr, every=20)
    for(i in 1:B) {
      pb(i)

      if(ngroup) {
        j <- integer(n)
        for(si in 1L : ngroup) {
          gi    <- group.inds[[si]]
          j[gi] <- sample(gi, length(gi), replace=TRUE)
        }
      }
      else j <- sample(1L : n, n, replace=TRUE)

      ## Note: If Strata is NULL, NULL[j] is still NULL
      f <- tryCatch(fitter(X[j,,drop=FALSE], Y[j,,drop=FALSE],
                           strata=Strata[j], ytarget=ytarget, ...),
                    error=function(...) list(fail=TRUE))
      if(length(f$fail) && f$fail) next

      cof <- f$coefficients
      if(any(is.na(cof))) next   # glm
      b <- b + 1L

      if(sc.pres) cof <- c(cof, 'log scale' = log(f$scale))
    
      non_repres_y <- length(vname) - length(cof)
      nry[i]       <- non_repres_y
      cof          <- process_ints(cof, num.intercepts(f), non_repres_y)
      
      if(coef.reps)             coefs[b, ] <- cof

      if(length(stat)) stats[b] <- f$stats[stat]

      bar <- bar + cof
      cof <- as.matrix(cof)
      cov <- cov + cof %*% t(cof)

      if(loglik) Loglik[b] <- oosl(f, matxv(X,cof), Y)
    }
  }
  else {
    clusterInfo <- list(name=deparse(substitute(cluster)))
    if(length(cluster) > n) {
      ## Missing obs were deleted during fit
      if(length(nac)) {
        j <- !is.na(naresid(nac, Y) %*% rep(1,ny))
        cluster <- cluster[j]
      }
    }

    if(length(cluster) != n)
      stop("length of cluster does not match # rows used in fit")

    if(any(is.na(cluster))) stop("cluster contains NAs")

    cluster <- as.character(cluster)

    clusters <- unique(cluster)
    nc <- length(clusters)
    Obsno <- split(1 : n, cluster)

    b <- 0
    pb <- setPb(B, type='Boot', onlytk=!pr, every=20)

    for(i in 1L : B) {
      pb(i)

      ## Begin addition Bill Pikounis
      if(ngroup) {
        j <- integer(0L)
        for(si in 1L : ngroup) {
          gi          <- group.inds[[si]]
          cluster.gi  <- cluster[gi]
          clusters.gi <- unique(cluster.gi)
          nc.gi       <- length(clusters.gi)
          Obsno.gci   <- split(gi, cluster.gi)
          j.gci       <- sample(clusters.gi, nc.gi, replace = TRUE)
          obs.gci     <- unlist(Obsno.gci[j.gci])
          j           <- c(j, obs.gci)
        }
        obs <- j
      }
      else {
        ## End addition Bill Pikounis (except for closing brace below)
        j   <- sample(clusters, nc, replace=TRUE)
        obs <- unlist(Obsno[j])
      }

      f <- tryCatch(fitter(X[obs,,drop=FALSE], Y[obs,,drop=FALSE],
                           strata=Strata[obs], ytarget=ytarget, ...),
                    error=function(...) list(fail=TRUE))
      if(length(f$fail) && f$fail) next

      cof <- f$coefficients

      if(any(is.na(cof))) next  # glm
      b <- b + 1L

      if(sc.pres) cof <- c(cof, 'log scale' = log(f$scale))

      non_repres_y <- length(vname) - length(cof)
      nry[i]       <- non_repres_y
      cof          <- process_ints(cof, num.intercepts(f), non_repres_y)

      if(coef.reps)    coefs[b,] <- cof
      if(length(stat)) stats[b] <- f$stats[stat]

      bar <- bar + cof
      cof <- as.matrix(cof)
      cov <- cov + cof %*% t(cof)
      if(loglik) Loglik[b] <- oosl(f, matxv(X,cof), Y)
    }
  }

  if(b < B) {
    warning('fit failure in ', B-b,
            ' resamples.  Consider specifying tol, maxit, opt_method, or other optimization criteria.')
    if(coef.reps) coefs <- coefs[1L : b,,drop=FALSE]
    Loglik <- Loglik[1L : b]
  }
  # if(nfit == 'orm') attr(coefs, 'intercepts') <- iref

  if(sum(nry) > 0) {
    cat('Counts of missing intercepts filled in by interpolation/extrapolation',
        '(median=', median(nry), 'out of', ns, 'intercepts)\n\n')
    print(table(nry))
  }
  
  bar           <- bar / b
  fit$B         <- b
  fit$seed      <- seed
  names(bar)    <- vname
  fit$boot.coef <- bar
  if(coef.reps) fit$boot.Coef <- coefs
  if(length(ytarget)) {
    fit$coefficients <- Cof
    fit$non.slopes   <- 1
    fit$interceptRef <- 1
    if(length(fit$linear.predictors) &&
       length(attr(fit$linear.predictors, 'intercepts')))
         attr(fit$linear.predictors, 'intercepts') <- 1
    fit$var <- infoMxop(fit$info.matrix, i=ikeep)
    for(i in 1 : length(fit$assign))
      fit$assign[[i]] <- fit$assign[[i]] - (ns - 1)
  }
  
  bar <- as.matrix(bar)
  cov <- (cov - b * bar %*% t(bar)) / (b - 1L)
  fit$orig.var <- fit$var
  # if(nfit == 'orm') attr(cov, 'intercepts') <- iref   1 if ytarget
  fit$var <- cov
  fit$info.matrix <- NULL
  fit$boot.loglik <- Loglik
  if(length(stat)) fit$boot.stats <- stats
  if(nfit == 'Rq') {
    newse <- sqrt(diag(cov))
    newt  <- fit$summary[, 1L] / newse
    newp  <- 2. * (1. - pt(abs(newt), fit$stats['n'] - fit$stats['p']))
    fit$summary[, 2L : 4L] <- cbind(newse, newt, newp)
  }
  if(length(clusterInfo)) clusterInfo$n <- nc
  fit$clusterInfo <- clusterInfo
  fit
}

bootplot <- function(obj, which=1 : ncol(Coef), X,
                     conf.int=c(.9,.95,.99),
                     what=c('density', 'qqnorm', 'box'),
                     fun=function(x) x,
                     labels., ...) {

  what <- match.arg(what)
  Coef <- obj$boot.Coef
  if(length(Coef) == 0) stop('did not specify "coef.reps=TRUE" to bootcov')

  Coef <- Coef[, which, drop=FALSE]

  if(! missing(X)) {
    if(! is.matrix(X)) X <- matrix(X, nrow=1)
    qoi <- matxv(X, Coef, bmat=TRUE)  # X %*% t(Coef)   ##nxp pxB = nxB
    if(missing(labels.)) {
      labels. <- dimnames(X)[[1]]
      if(length(labels.) == 0) {
        labels. <- as.character(1:nrow(X))
      }
    }
  } else {
    qoi <- t(Coef)
    nns <- num.intercepts(obj)
    if(missing(labels.)) {
      labels. <- paste(ifelse(which > nns, 'Coefficient of ', ''),
                       dimnames(Coef)[[2]], sep='')
    }
  }

  nq   <- nrow(qoi)
  qoi  <- fun(qoi)
  quan <- NULL

  if(what == 'box') {
    Co <- as.vector(Coef)
    predictor <- rep(colnames(Coef), each=nrow(Coef))
    p <- ggplot(data.frame(predictor, Co), aes(x=predictor, y=Co)) +
      xlab('Predictor') + ylab('Coefficient') +
      geom_boxplot() + facet_wrap(~ predictor, scales='free')
    return(p)
  }
  else if(what == 'density') {
    probs <- (1 + conf.int) / 2
    probs <- c(1 - probs, probs)
    quan <- matrix(NA, nrow=nq, ncol=2 * length(conf.int),
                   dimnames=list(labels., format(probs)))

      for(j in 1 : nq) {
        histdensity(qoi[j,], xlab=labels.[j], ...)
        quan[j,] <- quantile(qoi[j,], probs, na.rm=TRUE)
        abline(v=quan[j,], lty=2)
        title(sub=paste('Fraction of effects >', fun(0), ' = ',
                format(mean(qoi[j,] > fun(0))),sep=''), adj=0)
      }
  }
  else {
    for(j in 1 : nq) {
      qqnorm(qoi[j,], ylab=labels.[j])
      qqline(qoi[j,])
    }
  }

  invisible(list(qoi=drop(qoi), quantiles=drop(quan)))
}


## histdensity runs hist() and density(), using twice the number of
## class than the default for hist, and 1.5 times the width than the default
## for density

histdensity <- function(y, xlab, nclass, width, mult.width=1, ...) {
  y <- y[is.finite(y)]
  if(missing(xlab)) {
    xlab <- label(y)
    if(xlab == '') xlab <- as.character(sys.call())[-1]
  }

  if(missing(nclass)) nclass <- (logb(length(y),base=2)+1)*2

  hist(y, nclass=nclass, xlab=xlab, probability=TRUE, ...)
  if(missing(width)) {
    nbar <- logb(length(y), base = 2) + 1
    width <- diff(range(y))/nbar*.75*mult.width
  }

  lines(density(y,width=width,n=200))
  invisible()
}


confplot <- function(obj, X, against,
                     method=c('simultaneous', 'pointwise'),
                     conf.int=0.95,
                     fun=function(x) x,
                     add=FALSE, lty.conf=2, ...) {

  method <- match.arg(method)
  if(length(conf.int)>1) stop('may not specify more than one conf.int value')

  boot.Coef <- obj$boot.Coef
  if(length(boot.Coef) == 0) stop('did not specify "coef.reps=TRUE" to bootcov')

  if(!is.matrix(X)) X <- matrix(X, nrow=1)

  fitted <- fun(matxv(X, obj$coefficients))

  if(method == 'pointwise') {
    pred <- matxv(X, boot.Coef, bmat=TRUE)   ## n x B
    p <- fun(apply(pred, 1, quantile,
                   probs=c((1 - conf.int)/2, 1 - (1 - conf.int)/2),
                   na.rm=TRUE))
    lower <- p[1,]
    upper <- p[2,]
  }
  else {
    boot.Coef <- rbind(boot.Coef, obj$coefficients)
    loglik    <- obj$boot.loglik
    if(length(loglik) == 0) stop('did not specify "loglik=TRUE" to bootcov')

    crit  <- quantile(loglik, conf.int, na.rm=TRUE)
    qual  <- loglik <= crit
    boot.Coef <- boot.Coef[qual,,drop=FALSE]
    pred   <- matxv(X, boot.Coef, bmat=TRUE)  ## n x B
    upper  <- fun(apply(pred, 1, max))
    lower  <- fun(apply(pred, 1, min))
    pred   <- fun(pred)
  }

  if(!missing(against)) {
    lab <- label(against)
    if(lab == '') lab <- (as.character(sys.call())[-1])[3]

    if(add) lines(against, fitted, ...)
    else plot(against, fitted, xlab=lab, type='l', ...)

    lines(against, lower, lty=lty.conf)
    lines(against, upper, lty=lty.conf)
  }
  if(missing(against)) list(fitted=fitted, upper=upper, lower=lower)
  else invisible(list(fitted=fitted, upper=upper, lower=lower))
}

# Construct object suitable for boot:boot.ci
# Use boot package to get BCa confidence limits for a linear combination of
# model coefficients, e.g. bootcov results boot.Coef
# If boot.ci fails return only ordinary percentile CLs
bootBCa <- function(estimate, estimates, type=c('percentile','bca','basic'),
                    n, seed, conf.int=0.95) {
  type <- match.arg(type)
  if(type != 'percentile' && ! requireNamespace('boot', quietly = TRUE))
    stop('boot package not installed')
  estimate <- as.vector(estimate)
  ne <- length(estimate)
  if(!is.matrix(estimates)) estimates <- as.matrix(estimates)
  if(ncol(estimates) != ne)
    stop('no. columns in estimates != length of estimate')
  if(type == 'percentile') {
    a <- apply(estimates, 2, quantile,
                 probs=c((1-conf.int)/2, 1-(1-conf.int)/2), na.rm=TRUE)
    if(ne == 1) a <- as.vector(a)
    return(a)
  }
  lim <- matrix(NA, nrow=2, ncol=ne, dimnames=list(c('Lower','Upper'),NULL))
  R <- nrow(estimates)
  for(i in 1:ne) {
    w <- list(sim= 'ordinary',
              stype = 'i',
              t0 = estimate[i],
              t  = estimates[,i,drop=FALSE],
              R  = R,
              data    = 1:n,
              strata  = rep(1,   n),
              weights = rep(1/n, n),
              seed = seed,
              statistic = function(...) 1e10,
              call = match.call())

    cl <- try(boot::boot.ci(w, type=type, conf=conf.int), silent=TRUE)
    if(inherits(cl, 'try-error')) {
      cl <- c(NA,NA)
    warning('could not obtain bootstrap confidence interval')
    } else {
      cl <- if(type == 'bca') cl$bca else cl$basic
      m <- length(cl)
      cl <- cl[c(m - 1, m)]
    }
    lim[,i] <- cl
  }
  if(ne == 1) as.vector(lim) else lim
}
