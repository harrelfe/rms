lrm <- function(formula, data,subset, na.action=na.delete,
                method="lrm.fit", model=FALSE, x=FALSE, y=FALSE, 
                linear.predictors=TRUE, se.fit=FALSE, 
                penalty=0, penalty.matrix, tol=1e-7, strata.penalty=0,
                var.penalty=c('simple','sandwich'),
                weights, normwt=FALSE, scale=FALSE, ...)
{
  call <- match.call()
  var.penalty <- match.arg(var.penalty)
  m <- match.call(expand.dots=FALSE)
  mc <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action
  m$drop.unused.levels <- TRUE
  
  m[[1]] <- as.name("model.frame")
  nact <- NULL
  if(missing(data)) data <- NULL

  tform <- terms(formula, specials='strat', data=data)
  nstrata <- 1
  if(!missing(data) || (
						length(atl <- attr(tform,"term.labels")) && 
						any(atl!=".")))	{ ##X's present

    dul <- .Options$drop.unused.levels
    if(!length(dul) || dul) {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }

    X        <- Design(eval.parent(m))
    atrx     <- attributes(X)
    sformula <- atrx$sformula
    nact     <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms    <- atrx$terms
    attr(Terms, "formula") <- formula
    atr      <- atrx$Design
    mmcolnames <- atr$mmcolnames

    Y <- model.extract(X, 'response')
    offs <- atrx$offset
    if(!length(offs)) offs <- 0
    weights <- wt <- model.extract(X, 'weights')
    if(length(weights))
      warning('currently weights are ignored in model validation and bootstrapping lrm fits')
    if(model) m <- X
    stra <- attr(tform,'specials')$strat
    Strata <- NULL
    Terms.ns <- Terms
    if(length(stra)) {
      temp <- untangle.specials(Terms.ns, 'strat', 1)
      Terms.ns <- Terms.ns[-temp$terms]
      attr(Terms,   "factors") <- pmin(attr(Terms,"factors"),1)
      attr(Terms.ns,"factors") <- pmin(attr(Terms.ns,"factors"),1)
      Strata <- X[[stra]]
      nstrata <- length(levels(Strata))
    }
    X <- model.matrix(Terms.ns, X)
    alt <- attr(mmcolnames, 'alt')
    if(! all(mmcolnames %in% colnames(X)) && length(alt)) mmcolnames <- alt
    ## prn(colnames(X)); prn(mmcolnames)
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
      if(nrow(penalty.matrix)!=p || ncol(penalty.matrix)!=p) stop(
             paste("penalty.matrix does not have",p,"rows and columns"))
      psetup     <- Penalty.setup(atr, penalty)
      penalty    <- psetup$penalty
      multiplier <- psetup$multiplier
      if(length(multiplier)==1)
        penalty.matrix <- multiplier * penalty.matrix
      else {
        a <- diag(sqrt(multiplier))
        penalty.matrix <- a %*% penalty.matrix %*% a
      }
      }
  }
  else {
    X <- eval.parent(m)
    offs <- model.offset(X)
    if(! length(offs)) offs <- 0
    Y <- model.extract(X, 'response')
    Y <- Y[!is.na(Y)]
    Terms <- X <- NULL
    xpres <- FALSE
    penpres <- FALSE
    penalty.matrix <- NULL
    }  ##Model: y~. without data= -> no predictors
  
  if(method == "model.matrix") return(X)

  if(nstrata > 1) {
    if(scale) stop('scale=TRUE not implemented for stratified model')
    f <- lrm.fit.strat(X,Y,Strata,offset=offs,
                       penalty.matrix=penalty.matrix,
                       strata.penalty=strata.penalty,
                       tol=tol,
                       weights=weights,normwt=normwt, ...)
  }
  else {
    if(existsFunction(method)) {
      fitter <- getFunction(method)
      f <- fitter(X, Y, offset=offs,
                  penalty.matrix=penalty.matrix, tol=tol,
                  weights=weights, normwt=normwt, scale=scale, ...)
    }
    else stop(paste("unimplemented method:", method))
  }
  
  if(f$fail) stop("Unable to fit model using ", dQuote(method))
  
  f$call <- NULL
  if(model) f$model <- m
  if(x) {
    f$x <- X
    f$strata <- Strata
  }
  if(y) f$y <- Y
  nrp <- f$non.slopes
  if(penpres) {
    f$penalty <- penalty
    if(nstrata == 1) {
      ## Get improved covariance matrix
      v <- f$var
      if(var.penalty=='sandwich') f$var.from.info.matrix <- v
      f.nopenalty <- 
        fitter(X, Y, offset=offs, initial=f$coef, maxit=1, tol=tol,
               scale=scale)
      ##  info.matrix.unpenalized <- solvet(f.nopenalty$var, tol=tol)
      info.matrix.unpenalized <- f.nopenalty$info.matrix
      dag <- diag(info.matrix.unpenalized %*% v)
      f$effective.df.diagonal <- dag
      f$var <- if(var.penalty == 'simple') v else
      v %*% info.matrix.unpenalized %*% v
      df <- sum(dag[-(1:nrp)])
      lr <- f.nopenalty$stats["Model L.R."]
      pval <- 1 - pchisq(lr, df)
      f$stats[c('d.f.','Model L.R.','P')] <- c(df, lr, pval)  
    }
  }
  ass <- if(xpres) DesignAssign(atr, nrp, Terms) else list()
  
  if(xpres) {
    if(linear.predictors) names(f$linear.predictors) <- names(Y)
    else
        f$linear.predictors <- NULL

      if(se.fit) {
        if(nstrata > 1) stop('se.fit=T not implemented for strat')
        xint <- matrix(0, nrow=length(Y), ncol=f$non.slopes)
        xint[,1] <- 1
        X <- cbind(xint, X)
        se <- drop((((X %*% f$var) * X) %*% rep(1, ncol(X)))^.5)
        names(se) <- names(Y)
        f$se.fit <- se
      }
  }
  f <- c(f, list(call=call, Design=if(xpres)atr,
                 scale.pred=c("log odds","Odds Ratio"),
                 terms=Terms, assign=ass, na.action=nact, fail=FALSE,
                 interceptRef=1,
                 nstrata=nstrata, sformula=sformula))
  
  class(f) <- c("lrm","rms","glm")
  f
}

print.lrm <- function(x, digits=4, strata.coefs=FALSE, coefs=TRUE,
                      title='Logistic Regression Model', ...) {

  latex <- prType() == 'latex'
  
  z <- list()
  k <- 0
  
  if(length(x$freq) > 3) {
    k <- k + 1
    z[[k]] <- list(type='print', list(x$freq),
                   title='Frequencies of Responses')
  }
  if(length(x$sumwty)) {
    k <- k + 1
    z[[k]] <- list(type='print', list(x$sumwty),
                   title='Sum of Weights by Response Category')
  }
  if(!is.null(x$nmiss)) {  ## for backward compatibility
    k <- k + 1
    z[[k]] <- list(type='print', list(x$nmiss),
                   title='Frequencies of Missing Values Due to Each Variable')
  }
  else if(length(x$na.action)) {
    k <- k + 1
    z[[k]] <- list(type=paste('naprint',class(x$na.action),sep='.'),
                   list(x$na.action))
  }
  
  ns <- x$non.slopes
  nstrata <- x$nstrata
  if(!length(nstrata)) nstrata <- 1

  pm <- x$penalty.matrix
  penaltyFactor <- NULL
  if(length(pm)) {
    psc <- if(length(pm)==1) sqrt(pm)
    else
      sqrt(diag(pm))
    penalty.scale <- c(rep(0,ns),psc)
    cof <- matrix(x$coef[-(1:ns)], ncol=1)
    k <- k + 1
    z[[k]] <- list(type='print', list(as.data.frame(x$penalty, row.names='')),
                   title='Penalty factors')
    penaltyFactor <- as.vector(t(cof) %*% pm %*% cof)
  }

  ## ?ok to have uncommented next 3 lines?
  est.exp <- 1:ns
  if(length(x$est)) est.exp <-
    c(est.exp, ns + x$est[x$est + ns <= length(x$coefficients)])
  vv <- diag(x$var)
  cof <- x$coef
  if(strata.coefs) {
    cof <- c(cof, x$strata.coef)
    vv <- c(vv, vcov(x))
    ## TODO: implement in vcov:
    ## vv  <- c(vv,  vcov(x, which='strata.var.diag'))
    if(length(pm)) penalty.scale <- c(penalty.scale, rep(NA, x$nstrata - 1))
  }
  score.there <- nstrata==1 && (length(x$est) < length(x$coef) - ns)
  stats <- x$stats

  maxd <- stats['Max Deriv']
  ci <- x$clusterInfo
  misc <- reListclean(Obs   =stats['Obs'],
                      'Sum of weights'=stats['Sum of Weights'],
                      Strata=if(nstrata > 1) nstrata,
                      'Cluster on'  = ci$name,
                      'Clusters'    = ci$n,
                      'max |deriv|' = maxd)
  if(length(x$freq) < 4) {
    names(x$freq) <- paste(if(latex)'~~' else ' ',
                           names(x$freq), sep='')
    misc <- c(misc[1], x$freq, misc[-1])
  }
  lr   <- reListclean('LR chi2'    = stats['Model L.R.'],
                      'd.f.'       = round(stats['d.f.'], 3),
                      'Pr(> chi2)' = stats['P'],
                      Penalty      = penaltyFactor)
  disc <- reListclean(R2           = stats['R2'],
                      g            = stats['g'], gr=stats['gr'],
                      gp           = stats['gp'], Brier=stats['Brier'])
  discr <-reListclean(C       = stats['C'],
                      Dxy     = stats['Dxy'],
                      gamma   = stats['Gamma'],
                      'tau-a' = stats['Tau-a'])
  
  headings <- c('','Model Likelihood\nRatio Test',
                   'Discrimination\nIndexes',
                   'Rank Discrim.\nIndexes')
  
  data <- list(misc, c(lr, c(2,NA,-4,2)), c(disc,3), c(discr,3))
  k <- k + 1
  z[[k]] <- list(type='stats', list(headings=headings, data=data))

  if(coefs) {
    k <- k + 1
    z[[k]] <- list(type='coefmatrix',
                   list(coef=cof, se=sqrt(vv),
                        aux=if(length(pm)) penalty.scale,
                        auxname='Penalty Scale'))
  }
  
  if(score.there) {
    q <- (1:length(cof))[-est.exp]
    if(length(q)==1) vv <- x$var[q,q] else vv <- diag(x$var[q,q])
    Z <- x$u[q]/sqrt(vv)
    stats <- cbind(format(Z,digits=2), format(1-pchisq(Z^2,1),digits=4))
    dimnames(stats) <- list(names(cof[q]),c("Score Z","P"))
    k <- k + 1
    z[[k]] <- list(type='print', list(stats))
  }
  prModFit(x, title=title, z, digits=digits,
           coefs=coefs, ...)
}


