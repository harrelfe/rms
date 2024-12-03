lrm <- function(formula, data=environment(formula),
                subset, na.action=na.delete,
                method="lrm.fit", model=FALSE, x=FALSE, y=FALSE, 
                linear.predictors=TRUE, se.fit=FALSE, 
                penalty=0, penalty.matrix,
                var.penalty,
                weights, normwt=FALSE, scale, ...)
{
  call <- match.call()
  if(! missing(var.penalty)) warning('var.penalty is deprecated and ignored')
  if(! missing(scale))       warning('scale is deprecated and ignored; see lrm.fit transx=')

  callenv <- parent.frame()   # don't delay these evaluations
  weights <- if(! missing(weights)) eval(substitute(weights), data, callenv)
  subset  <- if(! missing(subset )) eval(substitute(subset),  data, callenv)

  data <-
    modelData(data, formula,
              weights=weights, subset=subset,
              na.action=na.action, callenv=callenv)

  tform <- terms(formula, data=data)
	if(length(atl <- attr(tform, "term.labels")) && 
						any(atl!="."))	{ ##X's present

    X        <- Design(data, formula=formula)

    atrx     <- attributes(X)
    nact     <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms    <- atrx$terms
    attr(Terms, "formula") <- formula
    sformula <- atrx$sformula
    atr      <- atrx$Design
    mmcolnames <- atr$mmcolnames

    Y <- model.extract(X, 'response')
    offs <- atrx$offset
    if(!length(offs)) offs <- 0
    weights <- wt <- model.extract(X, 'weights')
    if(length(weights))
      warning('currently weights are ignored in model validation and bootstrapping lrm fits')
    if(model) m <- X
    
    X <- model.matrix(Terms, X)
    alt <- attr(mmcolnames, 'alt')
    if(! all(mmcolnames %in% colnames(X)) && length(alt)) mmcolnames <- alt
    ## prn(colnames(X)); prn(mmcolnames)
    X <- X[, mmcolnames, drop=FALSE]

    colnames(X) <- atr$colnames
    xpres <- length(X) > 0

    p <- length(atr$colnames)
    n <- length(Y)

    penpres <- !(missing(penalty) && missing(penalty.matrix))
    
    if(!penpres) penalty.matrix <- matrix(0,ncol=p,nrow=p) else { 
      if(missing(penalty.matrix)) penalty.matrix <- Penalty.matrix(atr, X) else
      if(nrow(penalty.matrix)!=p || ncol(penalty.matrix)!=p) stop(
             paste("penalty.matrix does not have", p, "rows and columns"))
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
    X <- Design(data, formula=formula)

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

  if(! is.factor(Y)) Y <- as.vector(Y)  # in case Y is a matrix

  if(method != 'lrm.fit') stop('when not "model.frame" or "model.matrix" method must be "lrm.fit"')
    f <- lrm.fit(X, Y, offset=offs,
                 penalty.matrix=penalty.matrix,
                 weights=weights, normwt=normwt, ...)
  
  if(f$fail) {
    warning("Unable to fit model using ", dQuote(method))
    return(f)
    }
  
  f$call <- NULL
  if(model) f$model <- m
  if(x) f$x <- X
  if(y) f$y <- Y
  nrp <- f$non.slopes
  if(penpres) {
    f$penalty <- penalty
    ## Get improved covariance matrix
    v <- f$var
    # if(var.penalty=='sandwich') f$var.from.info.matrix <- v
    f.nopenalty <- 
        lrm.fit(X, Y, offset=offs, initial=f$coef, maxit=1, ...)
    ##  info.matrix.unpenalized <- solvet(f.nopenalty$var, tol=tol)
    info.matrix.unpenalized <- f.nopenalty$info.matrix
    dag <- diag(info.matrix.unpenalized %*% v)
    f$effective.df.diagonal <- dag
    f$var <- v
    # v %*% info.matrix.unpenalized %*% v
    df <- sum(dag[- (1 : nrp)])
    lr <- f.nopenalty$stats["Model L.R."]
    pval <- 1 - pchisq(lr, df)
    f$stats[c('d.f.','Model L.R.','P')] <- c(df, lr, pval)  
  }

  ass <- if(xpres) DesignAssign(atr, nrp, Terms) else list()
  
  if(xpres) {
    if(linear.predictors) names(f$linear.predictors) <- names(Y)
    else
        f$linear.predictors <- NULL

      if(se.fit) {
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
                 sformula=sformula,
                 terms=Terms, assign=ass, na.action=nact, fail=FALSE,
                 interceptRef=1))
  
  class(f) <- c("lrm", "rms", "glm")
  f
}

print.lrm <- function(x, digits=4, r2=c(0,2,4), coefs=TRUE,
                      pg=FALSE, title='Logistic Regression Model', ...) {

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
  stats <- x$stats

  maxd <- stats['Max Deriv']
  ci <- x$clusterInfo
  misc <- reListclean(Obs   =stats['Obs'],
                      'Sum of weights'=stats['Sum of Weights'],
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
                      Penalty      = penaltyFactor,
                      dec          = c(2, NA, -4, 2))
  newr2 <- grepl('R2\\(', names(stats))
  disc <- reListclean(R2        = if(0 %in% r2)  stats['R2'],
                      namesFrom = if(any(newr2)) stats[newr2][setdiff(r2, 0)],
                      g         = if(pg) stats['g'],
                      gr        = if(pg) stats['gr'],
                      gp        = if(pg) stats['gp'],
                      Brier     = stats['Brier'],
                      dec       = 3)

  discr <-reListclean(C       = stats['C'],
                      Dxy     = stats['Dxy'],
                      gamma   = stats['Gamma'],
                      'tau-a' = stats['Tau-a'],
                      dec     = 3)
  
  headings <- c('','Model Likelihood\nRatio Test',
                   'Discrimination\nIndexes',
                   'Rank Discrim.\nIndexes')
  
  data <- list(misc, lr, disc, discr)
  k <- k + 1
  z[[k]] <- list(type='stats', list(headings=headings, data=data))

  if(coefs) {
    k <- k + 1
    z[[k]] <- list(type='coefmatrix',
                   list(coef=cof, se=sqrt(vv),
                        aux=if(length(pm)) penalty.scale,
                        auxname='Penalty Scale'))
  }
  
  prModFit(x, title=title, z, digits=digits,
           coefs=coefs, ...)
}


