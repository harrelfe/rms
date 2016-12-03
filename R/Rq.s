## Modification of the rq function in the quantreg package written by
## Roger Koenker, Stephen Portnoy, Pin Tian Ng, Achim Zeileis,
## Philip Grosjean, Brian Ripley

Rq <- function (formula, tau = 0.5, data, subset, weights, na.action=na.delete, 
                method = "br", model = FALSE, contrasts = NULL,
                se='nid', hs=TRUE, x=FALSE, y=FALSE, ...) 
{
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), 
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  mf[[1]] <- as.name("model.frame")
  mf <- Design(eval.parent(mf))

  at <- attributes(mf)
  sformula <- at$sformula
  desatr <- at$Design
  attr(mf,'Design') <- NULL
  if (method == "model.frame") return(mf)
  mt <- at$terms
  weights <- model.weights(mf)
  Y <- model.response(mf)
  X <- model.matrix(mt, mf, contrasts)
  if (length(desatr$colnames)) colnames(X) <- c("Intercept", desatr$colnames)

  eps <- .Machine$double.eps^(2/3)
  Rho <- function(u, tau) u * (tau - (u < 0))
  if (length(tau) > 1)
    stop('does not allow more than one quantile to be estimated simultaneously')
  ## The following keeps quantreg from overriding latex generic in Hmisc
  ## library(quantreg, pos=length(search()) + 1)
  fit <- if (length(weights)) 
    quantreg::rq.wfit(X, Y, tau = tau, weights, method, ...)
  else quantreg::rq.fit(X, Y, tau = tau, method, ...)
  rownames(fit$residuals) <- rownames(dimnames(X)[[1]])
  rho <- sum(Rho(fit$residuals, tau))
  
  stats <- c(n=length(fit$residuals),
             p=length(fit$coefficients),
             g=GiniMd(fit$fitted.values),
             mad=mean(abs(fit$residuals), na.rm=TRUE))
  
  fit <- c(fit,
           list(
                na.action = at$na.action,
                formula   = formula,
                sformula  = sformula,
                terms     = mt,
                xlevels   = .getXlevels(mt, mf),
                call      = call,
                tau       = tau,
                method    = method,
                weights   = weights,
                residuals = drop(fit$residuals),
                rho       = rho,
#                fitted.values = drop(fit$fitted.values),
                model     = mf,
                Design    = desatr,
                assign    = DesignAssign(desatr, 1, mt),
                stats     = stats))
  attr(fit, "na.message") <- attr(m, "na.message")
  
  s <- quantreg::summary.rq(fit, covariance=TRUE, se=se, hs=hs)
  k <- s$coefficients
  nam <- names(fit$coefficients)
  rownames(k) <- nam
  fit$summary <- k
  cov <- s$cov
  dimnames(cov) <- list(nam, nam)
  fit$var <- cov
  fit$method <- method
  fit$se <- se
  fit$hs <- hs
  
  ## Remove the following since summary.rq has done its job
  if(!model) fit$model <- NULL
  if(!x) fit$x <- NULL else fit$x <- X[, -1, drop=FALSE]
  if(!y) fit$y <- NULL
  class(fit) <- c('Rq', 'rms',
                  if (method == "lasso") "lassorq"
                  else if (method == "scad") "scadrq",
                  "rq")
  fit
}

## Thanks to Duncan Murdoch for the formals alist substitute technique
RqFit <- function(fit, wallow=TRUE, passdots=FALSE)
  {
    w <- fit$weights
    if(length(w))
      {
        if(!wallow) stop('weights not implemented')
        g <- if(passdots) function(x, y, weights, tau, method, ...)
          quantreg::rq.wfit(cbind(Intercept=1., x), y, tau = tau, weights=weights,
                  method=method, ...)
        else function(x, y, weights, tau, method, ...)
          quantreg::rq.wfit(cbind(Intercept=1., x), y, tau = tau, weights=weights,
                  method=method)
        formals(g) <- eval(substitute(
                       alist(x=,y=, weights=,tau=deftau,method=defmethod,...=),
                       list(deftau=fit$tau, defmethod=fit$method)))
      }
    else
      {
        g <- if(passdots) function(x, y, tau, method, ...)
          quantreg::rq.fit(cbind(Intercept=1., x), y, tau = tau, method=method, ...)
        else
          function(x, y, tau, method, ...)
            quantreg::rq.fit(cbind(Intercept=1., x), y, tau = tau, method=method)
        formals(g) <-
          eval(substitute(alist(x=,y=, tau=deftau, method=defmethod,...=),
                          list(deftau=fit$tau, defmethod=fit$method)))
      }
    g
  }

print.Rq <- function(x, digits=4, coefs=TRUE, title, ...)
  {
    k <- 0
    z <- list()

    ftau <- format(round(x$tau, digits))
    if(missing(title))
      title <- if(prType() == 'latex')
        paste('Quantile Regression~~~~$\\tau$', ftau, sep='=') else
        paste('Quantile Regression\t\ttau:',     ftau)

    if(length(zz <- x$na.action))
      {
        k <- k + 1
        z[[k]] <- list(type=paste('naprint', class(zz)[1], sep='.'), list(zz))
      }
    
    s <- x$stats
    n <- s['n']; p <- s['p']; errordf <- n - p; g <- s['g']
    mad <- s['mad']

    ci <- x$clusterInfo
    misc <- reListclean(Obs=n, p=p, 'Residual d.f.'=errordf,
                     'Cluster on'=ci$name,
                     Clusters    =ci$n,
                     'mean |Y-Yhat|'=mad)
    disc <- reListclean(g=g)
    headings <- c('', 'Discrimination\nIndex')
    data     <- list(misc, c(disc,3))
    k <- k + 1
    z[[k]] <- list(type='stats', list(headings=headings, data=data))

    s <- x$summary
    k <- k + 1
    z[[k]] <- list(type='coefmatrix', 
                   list(coef = s[,'Value'],
                        se   = s[,'Std. Error'],
                        errordf = errordf))
    
    if (length(mes <- attr(x, "na.message")))
      {
        k <- k + 1
        z[[k]] <- list(type='cat', list(mes, '\n'))
      }

    prModFit(x, title=title, z, digits=digits, coefs=coefs, ...)
  }

latex.Rq <-
  function(object,
           file = paste(first.word(deparse(substitute(object))),
             ".tex", sep = ""), append=FALSE,
           which, varnames, columns=65, inline=FALSE, caption=NULL, ...)
{
  html <- prType() == 'html'
  if(html) file <- ''
  
  f   <- object
  tau <- f$tau
  at  <- f$Design
    
  w <- if (length(caption)) {
         if(html) paste('<div align=center><strong>', caption,
                      '</strong></div>', sep='')
         else
           paste("\\begin{center} \\bf", caption, "\\end{center}")
         }
    if (missing(which) & !inline)
      {
        Y <- paste("{\\rm ", as.character(formula(f))[2], 
                   "}", sep = "")
        w <- c(w, paste("\\[", Y, "_{", tau, "} = X\\beta, {\\rm \\ \\ where} \\\\ \\]", 
                        sep = ""))
      }
    if(missing(which)) which <- 1:length(at$name)
    if(missing(varnames)) varnames <- at$name
    
    cat(w, file = file, sep = if (length(w)) "\n"  else "", append = append)
    latexrms(f, file=file, append=TRUE, which=which, inline=inline,
             varnames=varnames, columns=columns, caption, ...)
  }

predict.Rq <- function(object, ..., kint=1, se.fit=FALSE)
  predictrms(object, ..., kint=kint, se.fit=se.fit)
