residuals.ols <-
  function(object, 
           type=c("ordinary","score","dfbeta","dfbetas","dffit","dffits","hat",
             "hscore","influence.measures"), ...)
{

  type <- match.arg(type)
  naa <- object$na.action

  if(type == 'influence.measures') {
    class(object) <- 'lm'
    return(influence.measures(object)$infmat)
    }
  
  if(type=="ordinary") return(naresid(naa, object$residuals))
  
  if(!length(object$x))stop("did not specify x=TRUE in fit")

  X <- cbind(Intercept=1, object$x)
  if(type=="score") return(naresid(naa, X * object$residuals))
  
  infl <- ols.influence(object)
  
  if(type=="hscore") return(naresid(naa, X *
       (object$residuals / (1 - infl$hat))))
  
  if(type=="dfbeta" | type=="dfbetas")
    {
      r <- t(coef(object) - t(coef(infl)))
      if(type=="dfbetas") r <- sweep(r, 2, diag(object$var)^.5, "/")
    }
  else
    if(type=="dffit") r <- (infl$hat * object$residuals)/(1 - infl$hat)
    else
      if(type=="dffits") r <- (infl$hat^.5)*object$residuals /
        (infl$sigma * (1 - infl$hat))
      else
        if(type=="hat") r <- infl$hat

  naresid(naa, r)
}

## lm.influence used to work but now it re-computes X for unknown
## reasons  24Nov00
ols.influence <- function(lm, x)
{
  GET <- function(x, what)
    {
      ## eventually, x[[what, exact=TRUE]]
      if(is.na(n <- match(what, names(x)))) NULL else x[[n]]
    }
  wt <- GET(lm, "weights")
  ## should really test for < 1/BIG if machine pars available
  e <- lm$residuals
  n <- length(e)
  if(length(wt))
    e <- e * sqrt(wt)
  beta <- lm$coef
  if(is.matrix(beta))
    {
      beta <- beta[, 1]
      e <- e[, 1]
      warning("multivariate response, only first y variable used")
    }
  na <- is.na(beta)
  beta <- beta[!na]
  p <- GET(lm, "rank")
  if(!length(p)) p <- sum(!na)
  R <- lm$qr$qr
  if(p < max(dim(R))) R <- R[1:p, 1:p]
  qr <- GET(lm, "qr")
  if(!length(qr)) {
    lm.x <- cbind(Intercept=1, GET(lm, "x"))
    if(length(wt))
      lm.x <- lm.x * sqrt(wt)
    if(any(na))
      lm.x <- lm.x[, !na, drop = FALSE]
    stop('not implemented')  # left.solve doesn't exist in R
    ## Q <- left.solve(R, lm.x)
  }
  else {
    if(length(wt) && any(zero <- wt == 0))
      {
        Q <- matrix(0., n, p)
        dimnames(Q) <- list(names(e), names(beta))
        Q[!zero,  ] <- qr.Q(qr)[, 1:p, drop = FALSE]
      }
    else
      {
        Q <- qr.Q(qr)
        if(p < ncol(Q))
          Q <- Q[, 1:p, drop = FALSE]
      }
  }
  h <- as.vector((Q^2 %*% array(1, c(p, 1))))
  h.res <- (1 - h)
  z <- e/h.res
  v1 <- e^2
  z <- t(Q * z)
  v.res <- sum(v1)
  v1 <- (v.res - v1/h.res)/(n - p - 1)
  ## BKW (2.8)
  dbeta <- backsolve(R, z)
  list(coefficients = t(beta - dbeta), sigma = sqrt(v1), hat = h)
}
