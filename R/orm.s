orm <- function(formula, data=environment(formula),
        subset, na.action=na.delete,
				method="orm.fit",
				family=c("logistic", "probit", "loglog", "cloglog", "cauchit"),
				model=FALSE, x=FALSE, y=FALSE,
				linear.predictors=TRUE, se.fit=FALSE,
				penalty=0, penalty.matrix,
        var.penalty=c('simple','sandwich'),
				scale=FALSE, maxit=30, weights, normwt=FALSE, ...)
{
  call        <- match.call()
  var.penalty <- match.arg(var.penalty)
  family      <- match.arg(family)

  nact        <- NULL

  tform <- terms(formula, data=data)
  if(! missing(data) || (
						length(atl <- attr(tform,"term.labels")) &&
						any(atl!=".")))	{ ##X's present

    callenv <- parent.frame()   # don't delay these evaluations
    subset  <- if(! missing(subset )) eval(substitute(subset),  data, callenv)
    weights <- if(! missing(weights)) eval(substitute(weights), data, callenv)

    X <-
      modelData(data, formula, weights=weights,
                subset = subset,
                na.action=na.action, callenv=callenv)

    X          <- Design(X, formula=formula)
    atrx       <- attributes(X)
    sformula   <- atrx$sformula
    nact       <- atrx$na.action
    if(method == "model.frame") return(X)
    Terms      <- atrx$terms
    attr(Terms, "formula") <- formula
    atr        <- atrx$Design
    mmcolnames <- atr$mmcolnames

    Y <- model.extract(X, 'response')
    offs <- atrx$offset
    if(!length(offs)) offs <- 0
    weights <- wt <- model.extract(X, 'weights')
    if(length(weights))
      warning('currently weights are ignored in model validation and bootstrapping orm fits')

    if(model) m <- X
    X <- model.matrix(sformula, X)
    alt <- attr(mmcolnames, 'alt')
    if(! all(mmcolnames %in% colnames(X)) && length(alt)) mmcolnames <- alt
    X <- X[, mmcolnames, drop=FALSE]
    colnames(X) <- atr$colnames

    xpres <- length(X) > 0

    p <- length(atr$colnames)
    n <- length(Y)

    penpres <- !(missing(penalty) && missing(penalty.matrix))
    if(penpres && missing(var.penalty))
      warning('default for var.penalty has changed to "simple"')

    if(!penpres) penalty.matrix <- matrix(0,ncol=p,nrow=p) else {
      if(missing(penalty.matrix)) penalty.matrix <- Penalty.matrix(atr, X) else
      if(nrow(penalty.matrix) != p || ncol(penalty.matrix) != p) stop(
             paste("penalty.matrix does not have",p,"rows and columns"))
      psetup <- Penalty.setup(atr, penalty)
      penalty <- psetup$penalty
      multiplier <- psetup$multiplier
      if(length(multiplier)==1)
        penalty.matrix <- multiplier*penalty.matrix
      else
        {
          a <- diag(sqrt(multiplier))
          penalty.matrix <- a %*% penalty.matrix %*% a
        }
    }
  }
  else
    {
      X <- eval.parent(m)
      offs <- model.offset(X)
      if(!length(offs)) offs <- 0
      Y <- model.extract(X, 'response')
      Y <- Y[!is.na(Y)]
      Terms <- X <- NULL
      xpres <- FALSE
      penpres <- FALSE
      penalty.matrix <- NULL
    }  ##Model: y~. without data= -> no predictors

  if(method=="model.matrix") return(X)

  if(existsFunction(method)) {
      fitter <- getFunction(method)
      f <- fitter(X, Y, family=family, offset=offs,
                  penalty.matrix=penalty.matrix,
                  scale=scale, maxit=maxit, weights=weights, normwt=normwt, ...)
    }
    else stop(paste("unimplemented method:", method))

  if(f$fail) {
    cat("Unable to fit model using ", dQuote(method), '\n')
    return(f)
  }

  f$call <- NULL
  f$sformula <- sformula
  if(model) f$model <- m
  if(x) f$x <- X
  if(y) f$y <- Y
  nrp  <- f$non.slopes
  info <- f$info.matrix
  if(penpres) {
    f$penalty <- penalty
    ## Get improved covariance matrix
    v <- infoMxop(info, invert=TRUE)
 
    if(var.penalty == 'sandwich') f$var.from.info.matrix <- v
    f.nopenalty <-
        fitter(X, Y, family=family, offset=offs, initial=f$coef, maxit=1,
               weights=weights, normwt=normwt)
    ##  info.matrix.unpenalized <- solvet(f.nopenalty$var, tol=tol)
    info.matrix.unpenalized <- infoMxop(f.nopenalty$info.matrix)
    ## Why can't just just come from f$info.matrix ??
    dag <- Matrix::diag(info.matrix.unpenalized %*% v)
    f$effective.df.diagonal <- dag
    f$var <- if(var.penalty == 'simple') v else
       v %*% info.matrix.unpenalized %*% v
    df   <- sum(dag[-(1:nrp)])
    lr   <- f.nopenalty$stats["Model L.R."]
    pval <- 1 - pchisq(lr, df)
    f$stats[c('d.f.','Model L.R.','P')] <- c(df, lr, pval)
  }
  ass <- if(xpres) DesignAssign(atr, nrp, Terms) else list()

  if(xpres) {
    if(linear.predictors) names(f$linear.predictors) <- names(Y)
    else
        f$linear.predictors <- NULL

      if(se.fit) {
        nx <- ncol(X)
        X  <- cbind(1, X)
        v  <- infoMxop(info, i=c(f$interceptRef, (nrp + 1) : (nrp + nx)), B=t(X))
        se <- drop(sqrt((t(v) * X) %*% rep(1, nx + 1)))
        names(se) <- names(Y)
        f$se.fit <- se
      }
  }
  f <- c(f, list(call=call, Design=if(xpres)atr,
                 scale.pred=if(f$family=='logistic')
                  c("log odds","Odds Ratio") else
                  if(f$family=='loglog') c("log hazard", "Hazard Ratio"),
                 terms=Terms, assign=ass, na.action=nact))
  class(f) <- c("orm","rms")
  f
}

print.orm <- function(x, digits=4, r2=c(0,2,4), coefs=TRUE, pg=FALSE,
                      intercepts=x$non.slopes < 10,
                      title, ...) {

  if(missing(title)) {
    title <- switch(x$family,
                    logistic = 'Logistic (Proportional Odds)',
                    probit = 'Probit',
                    cauchit = 'Cauchy',
                    loglog  = '-log-log',
                    cloglog = 'Complementary log-log')
    title <- paste(title, 'Ordinal Regression Model')
  }

  z <- list()
  k <- 0

  lf <- length(x$freq)
  if(lf > 3 && lf <= 20) {
    k <- k + 1
    z[[k]] <- list(type='print', list(x$freq),
                   title='Frequencies of Responses')
  }
  if(length(x$nmiss)) {  ## for backward compatibility
    k <- k + 1
    z[[k]] <- list(type='print', list(x$nmiss),
                   title='Frequencies of Missing Values Due to Each Variable')
  }
  else if(length(x$na.action)) {
    k <- k + 1
    z[[k]] <- list(type=paste('naprint', class(x$na.action), sep='.'),
                   list(x$na.action))
  }

  ns <- x$non.slopes
  ## coefficient intercepts kept: (fit.mult.impute)
  cik <- attr(x$coef, 'intercepts')  # esp. for fit.mult.impute
  if(length(cik) && intercepts) {
    warning('intercepts=TRUE not implemented for fit.mult.impute objects')
    intercepts <- FALSE
  }

  pm <- x$penalty.matrix
  penaltyFactor <- NULL
  if(length(pm)) {
    psc <- if(length(pm) == 1) sqrt(pm)
    else
      sqrt(diag(pm))
    penalty.scale <- c(rep(0, ns), psc)
    cof <- matrix(x$coef[-(1 : ns)], ncol=1)
    ## This logic does not handle fit.mult.impute objects
    k <- k + 1
    z[[k]] <- list(type='print', list(as.data.frame(x$penalty, row.names='')),
                   title='Penalty factors')
    penaltyFactor <- as.vector(t(cof) %*% pm %*% cof)
  }
  vv <- Matrix::diag(vcov(x, intercepts=if(intercepts) 'all' else 'none'))
  if(! intercepts) {
    nints <- if(!length(cik)) ns else {
      if(length(cik) == 1 && cik ==0) 0 else length(cik)
    }
    ints.to.delete <- if(ns == 0 || nints == 0) integer(0) else 1:nints
    vv <- c(rep(NA, nints), vv)
  }
  cof   <- x$coef
  stats <- x$stats

  maxd <- stats['Max Deriv']
  ci <- x$clusterInfo
  frq <- if(length(x$freq) < 4) {
    x$freq
    }

  misc <- reListclean(Obs           = stats['Obs'],
                      'Distinct Y'    = stats['Distinct Y'],
                      'Cluster on'  = ci$name,
                      Clusters      = ci$n,
                      'Median Y'    = stats['Median Y'],
                      'max |deriv|' = maxd)
  if(length(x$freq) < 4) {
    names(x$freq) <- paste(if(prType() == 'latex') '~~' else ' ',
                           names(x$freq), sep='')
    misc <- c(misc[1], x$freq, misc[-1])
  }
  lr   <- reListclean('LR chi2'    = stats['Model L.R.'],
                   'd.f.'       = round(stats['d.f.'],3),
                   'Pr(> chi2)' = stats['P'],
                   'Score chi2' = stats['Score'],
                   'Pr(> chi2)' = stats['Score P'],
                   Penalty      = penaltyFactor,
                   dec          = c(2,NA,-4,2,-4,2))
  newr2 <- grepl('R2\\(', names(stats))
  disc <- reListclean(R2=if(0 %in% r2) stats['R2'],
                      namesFrom = if(any(newr2)) stats[newr2][setdiff(r2, 0)],
                      g         = if(pg) stats['g'],
                      gr        = if(pg) stats['gr'],
                      '|Pr(Y>=median)-0.5|' = stats['pdm'],
                      dec       = 3)
  if(any(newr2)) names(disc)[names(disc) == 'R2m'] <- names(stats[newr2])

  discr <-reListclean(rho=stats['rho'], dec=3)

  headings <- c('',
                'Model Likelihood\nRatio Test',
                'Discrimination\n Indexes',
                'Rank Discrim.\nIndexes')
  data <- list(misc, lr, disc, discr)
  k <- k + 1
  z[[k]] <- list(type='stats', list(headings=headings, data=data))

  if(coefs) {
    k <- k + 1
    if(!intercepts) {
      j   <- - ints.to.delete
      cof <- cof[j]
      vv  <- vv[j]
      if(length(pm)) penalty.scale <- penalty.scale[j]
    }
    z[[k]] <- list(type='coefmatrix',
                   list(coef=cof, se=sqrt(vv),
                        aux=if(length(pm)) penalty.scale,
                        auxname='Penalty Scale'))
  }

  prModFit(x, title=title, z, digits=digits,
           coefs=coefs, ...)
}

Mean.orm <- function(object, codes=FALSE, ...)
  Mean.lrm(object, codes=codes, ...)

Quantile.orm <- function(object, codes=FALSE, ...)
{
  ns <- object$non.slopes
  if(ns < 2)
    stop('using this function only makes sense for >2 ordered response categories')
  if(codes) vals <- 1:length(object$freq)
  else {
    vals <- object$yunique
    if(!length(vals)) vals <- names(object$freq)
    vals <- as.numeric(vals)
    if(any(is.na(vals)))
      stop('values of response levels must be numeric for codes=FALSE')
  }
  f <- function(q=numeric(0), lp=numeric(0), X=numeric(0),
                intercepts=numeric(0), slopes=numeric(0),
                info=numeric(0), values=numeric(0),
                interceptRef=integer(0), trans=trans, conf.int=0,
                method=c('interpolated', 'discrete'))
  {
    inverse <- trans$inverse
    cumprob <- trans$cumprob
    deriv   <- trans$deriv
    ns <- length(intercepts)
    method <- match.arg(method)
    lp <- if(length(lp)) lp - intercepts[interceptRef] else matxv(X, slopes)
    lb <- matrix(sapply(intercepts, '+', lp), ncol = ns)
    if(method == 'interpolated'){
      m.yvals <- matrix(NA, nrow = nrow(lb), ncol = ns + 2)
      cp <- cbind(cumprob(lb), 0)
      for(j in 1:nrow(lb)){
        ws <- c(0, (cp[j, -ns-1] - cp[j, 1]) / (cp[j, ns] - cp[j, 1]), 1)
        m.yvals[j,] <- (1 - ws) * c(values[1], values) + ws * c(values, values[ns + 1])
      }
      z  <- sapply(1:nrow(lb),
                   function(i) approx(c(1, cp[i,]), m.yvals[i,],
                                      xout = cumprob(inverse(1 - q)), rule = 2)$y)
    }
    if(method == 'discrete'){
      m.cdf <- cbind(1 - cumprob(lb), 1)
      id <- apply(m.cdf, 1, FUN=function(x) {min(which(x >= q))[1]})
      z <- values[id]
    }
    names(z) <- names(lp)
    if(conf.int) {
      if(! length(X)) stop('when conf.int > 0 must provide X')
      lb.se <- matrix(NA, ncol = ns, nrow = nrow(X))
      # info.inverse <- infoMxop(info, invert=TRUE)
      idx <- which(names(c(intercepts, slopes)) %in% colnames(X))
      dlb.dtheta <- as.matrix(cbind(1, X))
      for(i in 1:ns){
        # v.i <- info.inverse[c(i, idx), c(i, idx)]
        # lb.se[, i] <- sqrt(diag(dlb.dtheta %*% v.i %*% t(dlb.dtheta)))
        # Compute (i, idx) portion of info inverse, multiplied by t(dlb.dtheta)
        v.i        <- infoMxop(info, i=c(i, idx), B=t(dlb.dtheta))
        lb.se[, i] <- sqrt(Matrix::diag(dlb.dtheta %*% v.i))
      }
      w <- qnorm((1 + conf.int) / 2)
      ci.ub <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] - w * lb.se[, i])}), ncol = ns)
      ci.lb <- matrix(sapply(1:ns, FUN=function(i) {1 - cumprob(lb[, i] + w * lb.se[, i])}), ncol = ns)
      if(method == 'interpolated'){
        z.ub <- sapply(1:nrow(lb),
                       function(i) approx(c(1, 1 - ci.lb[i,], 0), m.yvals[i,],
                                          xout = cumprob(inverse(1 - q)),
                                          rule = 2)$y)
        z.lb <- sapply(1:nrow(lb),
                       function(i) approx(c(1, 1 - ci.ub[i,], 0), m.yvals[i,],
                                          xout = cumprob(inverse(1 - q)),
                                          rule = 2)$y)
      }
      if(method == 'discrete'){
        id <- apply(cbind(ci.lb, 1), 1, FUN=function(x) {min(which(x >= q))[1]})
        z.ub <- values[id]
        id <- apply(cbind(ci.ub, 1), 1, FUN=function(x) {min(which(x >= q))[1]})
        z.lb <- values[id]
      }
      attr(z, 'limits') <- list(lower = z.lb,
                                upper = z.ub)
    }
    z
  }
  ## Re-write first derivative so that it doesn't need the f argument
  if(object$family == "logistic")
    object$trans$deriv <- function(x) {p <- plogis(x); p * (1. - p)}
  trans <- object$trans
  if(! length(trans)) trans <- probabilityFamilies$logistic
  ir <- object$interceptRef
  if(! length(ir)) ir <- 1
  formals(f) <- list(q=numeric(0), lp=numeric(0), X=numeric(0),
                     intercepts=object$coef[1:ns],
                     slopes=object$coef[-(1 : ns)],
                     info=object$info.matrix,
                     values=vals, interceptRef=ir, trans=trans,
                     conf.int=0, method=c('interpolated', 'discrete'))
  f
}




ExProb <- function(object, ...) UseMethod("ExProb")

ExProb.orm <- function(object, codes=FALSE, ...)
{
  ns <- object$non.slopes
  if(codes) vals <- 1:length(object$freq)
  else {
    vals <- as.numeric(object$yunique)
    if(any(is.na(vals)))
      stop('values of response levels must be numeric for codes=FALSE')
  }
  f <- function(lp=numeric(0), X=numeric(0), y=NULL,
                intercepts=numeric(0), slopes=numeric(0),
                info=numeric(0), values=numeric(0),
                interceptRef=integer(0), trans=trans, yname=NULL,
                conf.int=0)
  {
    cumprob <- trans$cumprob
    lp <- if(length(lp)) lp - intercepts[interceptRef] else matxv(X, slopes)
    prob <- cumprob(sapply(c(1e30, intercepts), '+', lp))
    dim(prob) <- c(length(lp), length(values))
    if(! length(y)) {
      colnames(prob) <- paste('Prob(Y>=', values, ')', sep='')
      y <- values
      result <- structure(list(y=values, prob=prob, yname=yname),
                       class='ExProb')
    }
    else {
      p <- t(apply(prob, 1,
                 function(probs) {
                   pa <- approx(values, probs, xout=y,
                                f=1, method='constant')$y
                   pa[y < min(values)] <- 1.
                   pa[y > max(values)] <- 0.
                   pa
                 } ) )
      if(length(y) > 1) {
        colnames(p) <- paste('Prob(Y>=', y, ')', sep='')
        p <- list(y=y, prob=p)
      }
      result <- structure(drop(p), class='ExProb')
    }
    if(conf.int){
      if(! length(X)) stop('must specify X if conf.int > 0')
      index <- sapply(y, FUN=function(x) {if(x <= min(values)) result <- 1
      else if(x >= max(values)) result <- length(values)
      else which(x <= values)[1] - 1})
      # info.inverse <- as.matrix(solve(info))
      idx <- which(names(c(intercepts, slopes)) %in% colnames(X))
      dlb.dtheta <- as.matrix(cbind(1, X))
      lb.se <- sapply(1:length(y), function(i)
        # diag(dlb.dtheta %*% info.inverse[c(index[i], idx), c(index[i], idx)] %*% t(dlb.dtheta))
        Matrix::diag(dlb.dtheta %*% infoMxop(info, i=c(index[i], idx), B=t(dlb.dtheta)))
        )
      lb.se <- matrix(sqrt(lb.se), ncol = length(y))
      m.alpha <- c(intercepts, slopes)[index]
      lb <- matrix(sapply(m.alpha, '+', lp), ncol = length(y))
      ci.ub <- matrix(sapply(1:length(y), FUN=function(i)
      {cumprob(lb[, i] + qnorm((1 + conf.int) / 2) * lb.se[, i])}), ncol = length(y))
      ci.lb <- matrix(sapply(1:length(y), FUN=function(i)
      {cumprob(lb[, i] - qnorm((1 + conf.int) / 2) * lb.se[, i])}), ncol = length(y))
      ci.ub[, which(y <= min(values))] <- ci.lb[, which(y <= min(values))] <- 1
      ci.ub[, which(y >= max(values))] <- ci.lb[, which(y >= max(values))] <- 0
      if(length(y) > 1)
        colnames(ci.ub) <- colnames(ci.lb) <- colnames(result$prob)
      attr(result, 'limits') <- list(lower = ci.lb,
                                     upper = ci.ub)
    }
    result
  }
  trans <- object$trans
  formals(f) <- list(lp=numeric(0), X=numeric(0), y=NULL,
                     intercepts=object$coef[1:ns],
                     slopes=object$coef[-(1 : ns)],
                     info=object$info.matrix, values=vals,
                     interceptRef=object$interceptRef,
                     trans=trans, yname=all.vars(object$terms)[1],
                     conf.int=0)

  f
}


plot.ExProb <- function(x, ..., data=NULL,
                        xlim=NULL, xlab=x$yname, ylab=expression(Prob(Y>=y)),
                        col=par('col'), col.vert='gray85',
                        pch=20, pch.data=21,
                        lwd=par('lwd'), lwd.data=lwd, lty.data=2, key=TRUE)
  {
    xlab <- xlab
    Y <- x[[2]]
    x <- x[[1]]
    if(!length(xlim)) xlim <- range(x)
    plot(0, 0, xlim=xlim, ylim=c(0,1), type='n', xlab=xlab, ylab=ylab, ...)
    if(!is.matrix(Y)) Y <- matrix(Y, nrow=1)
    xpts <- ypts <- numeric(0)
    for(i in 1:nrow(Y)) {
      y <- Y[i,]
      segments(x, y, x, c(y[-1], 0), col=col.vert)
      segments(x[-length(x)], y[-1], x[-1], y[-1], col=i, lwd=lwd)
      points(x, y, pch=pch, col=i)
      if(length(data)) { xpts <- c(xpts, x); ypts <- c(ypts, y) }
    }
    if(!length(data)) return(invisible())
    if(!is.list(data)) {
      groups <- rep(' ', length(data))
      Y      <- nomiss(data)
    } else {
      if(!is.list(data) || length(data) < 2) stop('inappropriate data=')
      data <- nomiss(data)
      groups <- data[[1]]
      Y      <- data[[2]]
    }
    i      <- 0
    for(g in unique(groups)) {
      i <- i + 1
      s <- groups == g
      y <- nomiss(Y[s])
      x <- sort(unique(y))
      f <- c(1., 1. - cumsum(table(y)) / length(y))
      if(x[1] > min(Y)) {
        x <- c(min(Y), x)
        f <- c(1.,     f)
      }
      y <- f[-length(f)]
      segments(x,             y,     x,     c(y[-1], 0),
               col=col.vert, lty=lty.data)
      segments(x[-length(x)], y[-1], x[-1], y[-1],
               col=i, lty=lty.data, lwd=lwd.data)
      points(x, y, pch=pch.data, col=i)
      xpts <- c(xpts, x); ypts <- c(ypts, y)
    }
    if(key && is.list(data))
      putKeyEmpty(xpts, ypts,
                  labels=unique(groups),
                  col=1:i, xlim=xlim, grid=FALSE)
    invisible()
  }
