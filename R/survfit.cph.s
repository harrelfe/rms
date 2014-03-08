survfit.cph <- function(formula, newdata, se.fit=TRUE, conf.int=.95, 
                        individual=FALSE, type=NULL, vartype=NULL,
                        conf.type=c('log', 'log-log', 'plain', 'none'),
                        id, ...) {
  object <- formula
  Call <- match.call()
  Call[[1]] <- as.name("survfit")  ## nicer output for the user
  censor <- FALSE

  type <- object$type
  if (! length(type)) {
    ## Use the appropriate one from the model
    w <- c("exact", "breslow", "efron")
    survtype <- match(object$method, w)
  }
  else {
    w <- c("kalbfleisch-prentice", "aalen", "efron",
           "kaplan-meier", "breslow", "fleming-harrington",
           "greenwood", "tsiatis", "exact")
    survtype <- match(match.arg(type, w), w)
    survtype <- c(1,2,3,1,2,3,1,2,3)[survtype]
  }
  vartype <- if(!length(vartype)) survtype
  else {
    w <- c("greenwood", "aalen", "efron", "tsiatis")
    vt <- match(match.arg(vartype, w), w)
    if(vt == 4) 2 else vt
  }
  
  if (!se.fit) conf.type <- "none"
  else conf.type <- match.arg(conf.type)

  xpres <- length(object$means) > 0
  y <- object$y
  if(!length(y)) stop('must use y=TRUE with fit')
  if(xpres) {
    X <- object$x
    if(!length(X)) stop('must use x=TRUE with fit')
    n <- nrow(X)
    xcenter <- object$means
    X <- X - rep(xcenter, rep.int(n, ncol(X)))
  }
  else {
    n <- nrow(y)
    X <- matrix(0, nrow=n, ncol=1)
  }

  strata <- object$Strata
  if(!length(strata)) strata <- rep(0,  n)
  offset <- object$offset
  if(!length(offset)) offset <- rep(0., n)
  weights <- object$weights
  if(!length(weights)) weights <- rep(1., n)

  missid <- missing(id)
  if (!missid) individual <- TRUE
  else if (missid && individual) id <- rep(0, n)
  else id <- NULL
  if (individual && attr(y, 'type') != "counting") 
    stop("The individual option is  only valid for start-stop data")
  
  ## Compute confidence limits for survival based on -log survival,
  ## constraining to be in [0,1]; d = std.error of cum hazard * z value
  ciupper <- function(surv, d) ifelse(surv==0, 0, pmin(1, surv*exp(d)))
  cilower <- function(surv, d) ifelse(surv==0, 0, surv*exp(-d))
  
  risk <- rep(exp(object$linear.predictors), length=n)
  ## need to center offset??
  ## coxph.fit centered offset inside linear predictors

  if(missing(newdata)) {
    X2 <- if(xpres) matrix(0., nrow=1, ncol=ncol(X)) else
    matrix(0., nrow=1, ncol=1)
    rq <- ro <- NULL
    newrisk <- 1
  }
  else {
    if (length(object$frail)) 
      stop("The newdata argument is not supported for sparse frailty terms")
    X2 <- predictrms(object, newdata, type='x', expand.na=FALSE)
    ## result with type='x' has attributes strata and offset which may be NULL
    rq <- attr(X2, 'strata')
    ro <- attr(X2, 'offset')
    n2 <- nrow(X2)
    if(length(rq) && any(levels(rq) %nin% levels(strata)))
      stop('new dataset has strata levels not found in the original')
    if(!length(rq)) rq <- rep(1,  n2)
    ro <- if(length(ro)) ro - mean(offset) else rep(0., n2)
    X2 <- X2 - rep(xcenter, rep.int(n2, ncol(X2)))
    newrisk <- exp(matxv(X2, object$coefficients) + ro)
  }

  y2 <- NULL
  if (individual) {
    if(missing(newdata))
      stop("The newdata argument must be present when individual=TRUE")
    isS <- sapply(newdata, is.Surv)
    if(sum(isS) != 1)
      stop("newdata must contain exactly one Surv object when individual=TRUE")
    y2 <- newdata[[which(isS)]]
  }
  g <- survfitcoxph.fit(y, X, weights, X2, risk, newrisk, strata,
                        se.fit, survtype, vartype,
                        if(length(object$var)) object$var else
                        matrix(0, nrow=1, ncol=1),
                        id=id, y2=y2, strata2=rq)
  if (!censor) {
    kfun <- function(x, keep) {
      if (is.matrix(x)) x[keep,, drop=FALSE] 
      else if (length(x) == length(keep)) x[keep] else x
    }
    keep <- g$n.event > 0
    if(length(g$strata)) {
      w <- factor(rep(names(g$strata), g$strata), names(g$strata))
      g$strata <- c(table(w[keep]))
    }
    g <- lapply(g, kfun, keep)
  }
  
  if (se.fit) {
    zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)
    if (conf.type=='plain') {
      u <- g$surv + zval* g$std.err * g$surv
      z <- g$surv - zval* g$std.err * g$surv
      g <- c(g, list(upper=pmin(u,1), lower=pmax(z,0),
                     conf.type='plain', conf.int=conf.int))
    }
    if (conf.type=='log')
      g <- c(g, list(upper=ciupper(g$surv, zval * g$std.err),
                     lower=cilower(g$surv, zval * g$std.err),
                     conf.type='log', conf.int=conf.int))
    if (conf.type=='log-log') {
      who <- (g$surv==0 | g$surv==1) #special cases
      xx <- ifelse(who, .1, g$surv)  #avoid some "log(0)" messages
      u <- exp(-exp(log(-log(xx)) + zval * g$std.err/log(xx)))
      u <- ifelse(who, g$surv + 0 * g$std.err, u)
      z <- exp(-exp(log(-log(xx)) - zval*g$std.err/log(xx)))
      z <- ifelse(who, g$surv + 0 * g$std.err, z)
      g <- c(g, list(upper=u, lower=z,
                     conf.type='log-log', conf.int=conf.int))
    }
  }
  g$requested.strata <- rq
  g$call <- Call
  class(g) <- c('survfit.cph', 'survfit')
  g
}
