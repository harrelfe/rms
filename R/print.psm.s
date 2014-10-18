print.psm <- function(x, correlation = FALSE, digits=4,
                      coefs=TRUE, latex=FALSE, title, ...)
{
  k <- 0
  z <- list()
  
  dist <- x$dist
  name <- survreg.distributions[[dist]]$name
  if(missing(title))
    title <- paste("Parametric Survival Model:", name, "Distribution")

  stats <- x$stats
  ci <- x$clusterInfo
  counts <- reVector(Obs   = stats['Obs'],
                     Events= stats['Events'],
                     'Cluster on' = ci$name,
                     Clusters = ci$n,
                     'Sum of Weights'=stats['Sum of Weights'],
                     sigma = if(length(x$scale) == 1) x$scale)
  lr <- reVector('LR chi2'     = stats['Model L.R.'],
                 'd.f.'        = stats['d.f.'],
                 'Pr(> chi2)'  = stats['P'])
  disc <- reVector(R2=stats['R2'], Dxy=stats['Dxy'],
                   g=stats['g'],   gr=stats['gr'])

  headings <- list('',
                   c('Model Likelihood','Ratio Test'),
                   c('Discrimination','Indexes'))
  digcounts <- c(NA, NA, NA,
                 if(length(ci$name)) NA,
                 if(length(ci$n))    NA,
                 if(length(x$scale) == 1) 4)

  data <- list(c(counts, digcounts), c(lr, c(2,NA,-4)), c(disc, 3))
  k <- k + 1
  z[[k]] <- list(type='stats', list(headings=headings, data=data))

  summary.survreg <- getS3method('summary', 'survreg')
  if(!x$fail) x$fail <- NULL    # summary.survreg uses NULL for OK
  s <- summary.survreg(x, correlation=correlation)
  k <- k + 1
  z[[k]] <- list(type='coefmatrix',
                 list(coef = s$table[,'Value'],
                      se   = s$table[,'Std. Error']))

  if (correlation && length(correl <- s$correlation)) {
    p <- ncol(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits = digits))
      correl[!ll] <- ""
      k <- k + 1
      z[[k]] <- list(type='print',
                     list(correl[-1, -p, drop = FALSE], quote = FALSE),
                     title='Correlation of Coefficients')
    }
  }
  
  prModFit(x, title=title, z, digits=digits, coefs=coefs, latex=latex, ...)
  invisible()
}

#  wt <- x$weights
#  fparms <- x$fixed
#  coef <- c(x$coef, x$parms[!fparms])
#  resid <- x$residuals
#  dresid <- x$dresiduals
#  n <- length(resid)
#  p <- x$rank
#  if(!length(p)) p <- sum(!is.na(coef))
#  if(!p)
#    {
#      warning("This model has zero rank --- no summary is provided")
#      return(x)
#    }
#  nsingular <- length(coef) - p
#  rdf <- x$df.resid
#  if(!length(rdf))
#    rdf <- n - p
#  R <- x$R   #check for rank deficiencies
#  if(p < max(dim(R)))
#    R <- R[1:p,     #coded by pivoting
#           1:p]
#  if(length(wt))
#    {
#      wt <- wt^0.5
#      resid <- resid * wt
#      excl <- wt == 0
#      if(any(excl))
#        {
#          warning(paste(sum(excl), 
#                        "rows with zero weights not counted"))
#          resid <- resid[!excl]
#          if(!length(x$df.residual))
#            rdf <- rdf - sum(excl)
#        }
#    }
#  famname <- x$family["name"]
#  if(!length(famname)) famname <- "Gaussian"
#  scale <- x$fparms
#  nas <- is.na(coef)
#  cnames <- names(coef[!nas])
#  coef <- matrix(rep(coef[!nas], 4), ncol = 4)
#  dimnames(coef) <- list(cnames, c("Value", "Std. Error", "z value", "p"))
#  stds <- sqrt(diag(x$var[!nas,!nas,drop=FALSE]))
#  coef[, 2] <- stds
#  coef[, 3] <- coef[, 1]/stds
#  coef[, 4] <- 2*pnorm(-abs(coef[,3]))
#  if(correlation)
#    {
#      if(sum(nas)==1) ss <- 1/stds else ss <- diag(1/stds)
#      correl <- ss %*% x$var[!nas, !nas, drop=FALSE] %*% ss
#      dimnames(correl) <- list(cnames, cnames)
#    }
#  else
#    correl <- NULL
#  ocall <- x$call
#  if(length(form <- x$formula))
#    {
#      if(!length(ocall$formula))
#        ocall <- match.call(get("survreg"), ocall)
#      ocall$formula <- form
#    }
#  dig <- .Options$digits
#  survival:::print.summary.survreg(
#                        list(call = ocall, terms = x$terms, coefficients = coef#,
#                             df = c(p, rdf), deviance.resid = dresid,
#                             var=x$var, correlation = correl, deviance = devian#ce(x),
#                             null.deviance = x$null.deviance, loglik=x$loglik,
#                             iter = x$iter,
#                             nas = nas))
#  options(digits=dig)   #recovers from bug in print.summary.survreg
#  invisible()
#}

## Mod of print.summary.survreg from survival5 - suppresses printing a
## few things, added correlation arg

print.summary.survreg2 <-
  function (x, digits = max(options()$digits - 4, 3),
            correlation=FALSE, ...) 
  {
    correl <- x$correl
    n <- x$n
    if (is.null(digits)) 
      digits <- options()$digits
    print(x$table, digits = digits)
    if (nrow(x$var) == length(x$coefficients)) 
      cat("\nScale fixed at", format(x$scale, digits = digits), 
          "\n")
    else
      if (length(x$scale) == 1) 
        cat("\nScale=", format(x$scale, digits = digits), "\n")
      else {
        cat("\nScale:\n")
        print(x$scale, digits = digits, ...)
      }
    
    if (correlation && length(correl)) {
      p <- dim(correl)[2]
      if (p > 1) {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits = digits))
        correl[!ll] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
    cat("\n")
    invisible(NULL)
  }
