orm.fit <- function(x=NULL, y,
                    family='logistic',
                    offset=0., initial, 
                    maxit=12L, eps=.005, tol=1e-7, trace=FALSE,
                    penalty.matrix=NULL, scale=FALSE)
{	
  cal <- match.call()
  len.penmat <- length(penalty.matrix)

  ## Extreme value type I dist = Gumbel maximum = exp(-exp(-x)) = MASS:::pgumbel
  ## Gumbel minimum = 1 - exp(-exp(x))
  families <-
    list(logistic =
         list(cumprob=function(x)    1 / (1 + exp(-x)),
              inverse=function(x)    log(x / (1 - x)),
              deriv  =function(x, f) f * (1 - f),
              deriv2 =function(x, f, deriv) f * (1 - 3*f + 2*f*f) ),
         probit =
         list(cumprob=pnorm,
              inverse=qnorm,
              deriv  =function(x, ...)      dnorm(x),
              deriv2 =function(x, f, deriv) - deriv * x),
         loglog =
         list(cumprob=function(x)      exp(-exp(-x)),
              inverse=function(x)      -log(-log(x    )),
              deriv  =function(x, ...) exp(-x - exp(-x)),
              deriv2 =function(x, ...) ifelse(abs(x) > 200, 0,
                exp(-x - exp(-x)) * (-1 + exp(-x)))),
         cloglog =
         list(cumprob=function(x)      1 - exp(-exp(x)),
              inverse=function(x)      log(-log(1 - x)),
              deriv  =function(x, ...) exp( x - exp( x)),
              deriv2 =function(x, f, deriv) ifelse(abs(x) > 200, 0,
                deriv * ( 1 - exp( x)))),
         cauchit =
         list(cumprob=pcauchy, inverse=qcauchy,
              deriv  =function(x, ...) dcauchy(x),
              deriv2 =function(x, ...) -2 * x * ((1 + x*x)^(-2)) / pi))
                     
  ## Check:
  ## P(x) = plogis(x); P'(x) = P(x) - P(x)^2
  ## d <- function(x) plogis(x) - 3*plogis(x)^2 + 2*plogis(x)^3
  ## x <- seq(-3, 3, length=150)
  ## plot(x, d(x), type='l')
  ## ad <- c(NA,diff(dlogis(x))/(x[2]-x[1]))
  ## lines(x, ad, col='red')

  familiesDefined <- names(families)
  sfam  <- substitute(family)
  csfam <- as.character(sfam)
  if(length(csfam) == 1 && csfam %in% familiesDefined) {
    fam <- families[[csfam]]
    family <- csfam
  }
  else if(is.character(family) && family %in% familiesDefined)
    fam <- families[[family]]
  else {
    fam <- family
    family <- fam$name
  }
  
  n <- length(y)
  
  initial.there <- !missing(initial)
  if(! length(x)) {
      nx <- 0
      xname <- NULL
      x <- 0
    }
  else  {
      if(! is.matrix(x)) x <- as.matrix(x)
      dx <- dim(x)
      nx <- dx[2L]
      if(dx[1] != n) stop("x and y must have same length")
      xname <- dimnames(x)[[2]]
      if(!length(xname)) xname <- paste("x[", 1 : nx, "]", sep="")
      if(scale) {
        x <- scale(x)
        scinfo <- attributes(x)[c('scaled:center', 'scaled:scale')]
        xbar <- scinfo[[1]]
        xsd  <- scinfo[[2]]
      }
    }


  ynumeric <- is.numeric(y)
  if(ynumeric) {
    mediany <- quantile(y, probs=.5, type=1L)
    yu <- sort(unique(y))
    kmid <- max(1, which(yu == mediany) - 1L)
  }
  # For large n, as.factor is slow
  # if(!is.factor(y)) y <- as.factor(y)
  if(is.factor(y)) {
    ylevels <- levels(y)
    y       <- unclass(y)
  }
  else {
    ylevels <- sort(unique(y))
    y       <- match(y, ylevels)
  }
  if(!ynumeric) {
    mediany <- quantile(y, probs=.5, type=1L)
    kmid    <- max(1, which(1L : length(ylevels) == mediany) - 1L)
  }

  kint <- length(ylevels) - 1L
  if(kint == 1) kmid <- 1
  ofpres <- !all(offset == 0)
  if(ofpres && length(offset) != n) stop("offset and y must have same length")
  
  if(n < 3) stop("must have >=3 non-missing observations")
  numy        <- tabulate(y)
  names(numy) <- ylevels
  p           <- as.integer(nx + kint)
  
  if(missing(initial)) {
      cp   <- (n - cumsum(numy)[- length(numy)]) / n
      names(cp) <- NULL
      initial <- fam$inverse(cp)
      if(ofpres) initial <- initial - mean(offset)
  }
  if(length(initial) < p)
      initial <- c(initial, rep(0, p - length(initial)))
  
  loglik <- -2 * sum(numy * log(numy / n))

  if(len.penmat) {
    if(nx > 0) {
      if(len.penmat == 0) penalty.matrix <- matrix(0, nrow=nx, ncol=nx)
      if(nrow(penalty.matrix) != nx || ncol(penalty.matrix) != nx) 
        stop(paste("penalty.matrix does not have", nx, "rows and columns"))
      penmat <- rbind(
        matrix(0, ncol=kint+nx, nrow=kint),
        cbind(matrix(0, ncol=kint, nrow=nx), penalty.matrix))
    }
    else penmat <- matrix(0, ncol=kint, nrow=kint)
  }
  else penmat <- NULL
  
  if(nx==0 & !ofpres) {
      loglik <- rep(loglik, 2)
      z <- list(coef=initial, u=rep(0,kint))
  }
  if(ofpres) {
      ##Fit model with only intercept(s) and offset
      z <- ormfit(NULL, y, kint, 0, initial, offset=offset,
                  maxit=maxit, tol=tol, eps=eps, trace=trace, fam)
      if(z$fail) return(structure(list(fail=TRUE), class="orm"))
      loglik <- c(loglik, z$loglik)
      initial <- z$coef
  }
 
  if(nx > 0) {
      ##Fit model with intercept(s), offset, covariables
      z <- ormfit(x, y, kint, nx, initial=initial, offset=offset, penmat=penmat,
                  maxit=maxit, tol=tol, eps=eps, trace=trace, fam)
      if(z$fail) return(structure(list(fail=TRUE), class="orm"))
      loglik <- c(loglik, z$loglik)
      kof  <- z$coef
      ## Compute linear predictor before unscaling beta, as x is scaled
      lp <- matxv(x, kof, kint=kmid)
      
      info <- z$v
      if(scale) {
        attr(info, 'scale') <- list(mean=xbar, sd=xsd)
        betas <- kof[- (1 : kint)]
        kof[1 : kint] <- kof[1 : kint] - sum(betas * xbar / xsd)
        kof[-(1 : kint)] <- betas / xsd
      }
  } else lp <- rep(kof[kmid], n)


  ## Keep variance matrix for middle intercept and all predictors
  ## Middle intercept take to be intercept corresponding to y that is
  ## closest to the median y

  i <- if(nx > 0) c(kmid, (kint + 1):p) else kmid
  v <- tryCatch(as.matrix(solve(info, tol=tol)[i, i]))
  if(inherits(v, 'try-error')) {
    cat('Singular information matrix\n')
    return(structure(list(fail=TRUE), class="orm"))
  }
  if(scale) {
    trans <- rbind(cbind(1, matrix(0, nrow=1, ncol=nx)),
                   cbind(-matrix(rep(xbar/xsd, 1), ncol=1), diag(1 / xsd)))
    v <- t(trans) %*% v %*% trans
  }
  name <- if(kint == 1) "Intercept" else
    paste("y>=", ylevels[-1L], sep="")
  name       <- c(name, xname)
  names(kof) <- name

  dimnames(v) <- list(name[i], name[i])
  if(kint > 1L) attr(v, 'intercepts') <- kmid

  llnull   <- loglik[length(loglik) - 1L]
  model.lr <- llnull - loglik[length(loglik)]
  model.df <- nx
  if(initial.there)
    model.p <- score <- score.p <- NA
  else {
    score <- z$score
    if(model.df == 0)
      model.p <- score.p <- 1.
    else {
      model.p <- 1. - pchisq(model.lr, model.df)
      score.p <- 1. - pchisq(score,    model.df)
    }
  }
  
  r2     <- 1. - exp(-model.lr / n)
  r2.max <- 1. - exp(-llnull / n)
  r2     <- r2 / r2.max
  if(kint > 1L) attr(lp, 'intercepts') <- kmid
  g  <- GiniMd(lp)
  ## compute average |difference| between 0.5 and the condition
  ## probability of being >= marginal median
  pdm <- mean(abs(fam$cumprob(lp) - 0.5))
  rho <- cor(rank(lp), rank(y))
  ## Somewhat faster: 
  ## rho <- .Fortran('rcorr', cbind(lp, y), as.integer(n), 2L, 2L, r=double(4),
  ##                 integer(4), double(n), double(n), double(n), double(n),
  ##                 double(n), integer(n), PACKAGE='Hmisc')$r[2]
  
  stats <- c(n, length(numy), mediany, z$dmax, model.lr, model.df,
             model.p, score, score.p, rho, r2, g, exp(g), pdm)
  
  nam <- c("Obs", "Unique Y", "Median Y", "Max Deriv",
           "Model L.R.", "d.f.", "P", "Score", "Score P",
           "rho", "R2", "g", "gr", "pdm") 
  names(stats) <- nam
  
  retlist <- list(call=cal, freq=numy, yunique=ylevels,
                  stats=stats, fail=FALSE, coefficients=kof,
                  var=v, ## u=z$u,
                  family=family, trans=fam,
                  deviance=loglik,
                  non.slopes=kint,
                  interceptRef=kmid,
                  linear.predictors=lp,
                  penalty.matrix=if(nx > 0 && any(penalty.matrix != 0))
                  penalty.matrix else NULL,
                  info.matrix=info)
  
  class(retlist) <- 'orm'
  retlist
}

ormfit <- function(x, y, kint, nx, initial, offset, penmat=NULL,
                   maxit=12L, eps=.005, tol=1e-7, trace=FALSE, fam) {
  n <- length(y)
  p <- as.integer(kint + nx)
  ymax <- kint + 1L
  iter <- 0L
  oldL <- 1e100
  coef <- initial
  del  <- rep(0., p)
  curstp <- 1.

  f   <- fam$cumprob
  fp  <- fam$deriv
  fpp <- fam$deriv2
  
  ep <- .Machine$double.eps
  fptest  <- function(x) fp (x, f(x))
  fpptest <- function(x) fpp(x, f(x), fp(x, f(x)))
  
  repeat {
    if(iter >= maxit) {
      cat('Did not converge in', maxit, 'iterations\n')
      return(list(fail=TRUE))
    }
    iter <- iter + 1L
    ## Compute linear predictor less intercept
    xb <- offset + (if(nx == 0L) 0. else x %*% coef[-(1L : kint)])
    ## Compute current Prob y=observed y
    ## P <- rep(0., n)
    ## P[y == 1]           <- f(xb + coef[1]) - f(xb + coef[2])
    ## if(ymax > 2) P[y > 1 & y < ymax] <- f(xb + coef[y]) - f(xb + coef[y+1])
    ## P[y == ymax]        <- f(xb + coef[kint])
    ints <- c(1e100, coef[1:kint], -1e100)
    xby <- xb + ints[y]; xby1 <- xb + ints[y + 1L]
    fa <- f(xby)
    fb <- f(xby1)
    P  <- fa - fb
    
    ## Compute -2 log likelihood
    L <- -2. * sum(log(P))
    ## Compute components of 1st and 2nd derivatives
    fpa  <- fp(xby,   fa)
    fpb  <- fp(xby1,  fb)
    fppa <- fpp(xby,  fa, fpa)
    fppb <- fpp(xby1, fb, fpb)
    
    if(abs(fptest(ints[1]      + max(xb))) > ep ||
       abs(fptest(ints[kint+2] + min(xb))) > ep)  stop('logic error 1')
    if(abs(fpptest(ints[1]      + max(xb))) > ep ||
       abs(fpptest(ints[kint+2] + min(xb))) > ep) stop('logic error 2')

    slow <- FALSE
    
    ## To compute score vector and info matrix in R instead of Ratfor/Fortran
    if(slow) {
      U <- rep(0., p)
      for(m in 1:kint)
        for(j in 1:n) U[m] <- U[m] +
          (fpa[j]*(y[j]-1==m) - fpb[j]*(y[j]==m)) / P[j]
      if(nx > 0) for(m in (kint+1):p)
        for(j in 1:n) U[m] <- U[m] +
          (fpa[j] - fpb[j]) * x[j,m-kint] / P[j]
      
      ## To compute the (wastefully sparse) information matrix in R:
      V <- matrix(NA, p, p)
      for(m in 1:p) {
        for(k in 1:m) {
          V[m,k] <- 0
          for(j in 1:n) {
            z <- y[j]
            pa  <- fpa[j];  pb  <- fpb[j]
            ppa <- fppa[j]; ppb <- fppb[j]
            w <- 1/(P[j]*P[j])
            v <- NA  ## temp
            v <-
              if(m <= kint & k <= kint)
                - w * (pa*(z-1==m) - pb*(z==m)) *
                  (pa*(z-1==k) - pb*(z==k)) +
                    (ppa*(z-1==m)*(m==k) -
                     ppb*(z==m)*(m==k))/P[j]
              else if(m > kint & k <= kint)
                x[j,m-kint] / P[j] *
                  (-1/P[j] * (pa - pb) * (pa*(z-1==k) - pb*(z==k)) +
                   ppa*(z-1==k) - ppb*(z==k))
              else if(m > kint & k > kint)
                x[j,m-kint] * x[j,k-kint] / P[j] *
                  (-1/P[j] * (pa - pb) * (pa-pb) + ppa - ppb)
              else {cat('nd\n');prn(c(m,k));9999}
            V[m,k] <- V[m,k] + v
          }}}
      V[col(V) > row(V)] <- t(V)[row(V) < col(V)]
      V <- -V
    } else {
      l <- if(kint == 1) p ^ 2
      else
        as.integer(nx * nx + 2L * kint * nx + 3L * kint - 2L)
      lia <- as.integer(p + 1L)
      w <- .Fortran('ormuv', n, p, kint, nx, x, y, P, fpa, fpb, fppa, fppb,
                    u=double(p), v=double(l), ja=integer(l), ia=integer(lia),
                    l=as.integer(l), lia=lia, integer(p), PACKAGE='rms')
      U <-  w$u
      V <- if(any(is.nan(w$v))) NULL else {  ## assume step-halving happens
        V <- if(kint == 1L) matrix(w$v, nrow=p, ncol=p)
        else new('matrix.csr', ra=w$v, ja=w$ja, ia=w$ia, dimension=c(p, p))
        V <- (V + t(V))/2.   # force symmetry; chol() complains if 1e-15 off
      }
    }
    dmax <- max(abs(U))
    ## Don't try to step halve if initial estimates happen to be
    ## MLEs; allow for a tiny worsening of LL without step-halving if
    ## max absolute first derivative is small
    if(trace) cat('-2logL=', L, ' step=', curstp, ' delta LL=', oldL-L,
                  ' max |deriv|=', dmax, '\n')
    if(dmax < 1e-9 && abs(L - oldL) < 0.01*eps) break
    if(abs(oldL - L) < eps) break
    if(L > oldL) {  # going the wrong direction
      if(iter == 1L) {
        cat('First iteration moved parameter estimates in wrong direction.\n')
        return(list(fail=TRUE))
      }
      ## Try step-halving
      curstp <- curstp / 2.
      coef <- coef - curstp * del
      next
    }
    curstp <- 1.
    del <- tryCatch(solve(V, U, tol=tol))
    if(inherits(del, 'try-error')) {
      cat('Singular information matrix\n')
      return(list(fail=TRUE))
    }
    if(iter == 1L) score <- U %*% del
    coef <- coef + del
    oldL <- L
  }
  ## Converged
  list(coef=coef, v=V, loglik=L, score=score, dmax=dmax,
       iter=iter, fail=FALSE, eps=eps)
}
