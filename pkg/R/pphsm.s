pphsm <- function(fit)
{
  warning("at present, pphsm does not return the correct covariance matrix")

  clas <- c(oldClass(fit), fit$fitFunction)
  if(!any(c('psm','survreg') %in% clas))
    stop("fit must be created by psm or survreg")
  if(fit$dist %nin% c('exponential','weibull'))
    stop("fit must have used dist='weibull' or 'exponential'")

  fit$coefficients <- -fit$coefficients/fit$scale
  fit$scale.pred   <- c("log Relative Hazard","Hazard Ratio")
  oldClass(fit) <- c("pphsm",oldClass(fit))

  fit
}

print.pphsm <- function(x, digits = max(options()$digits - 4, 3),
                        correlation = TRUE, ...)
{
  if (length(f <- x$fail) && f)
      stop(" Survreg failed.  No summary provided")
  
  cat("Parametric Survival Model Converted to PH Form\n\n")

  stats <- x$stats
  stats[3] <- round(stats[3],2)
  stats[5] <- round(stats[5],4)
  stats[6] <- round(stats[6],2)
  print(format.sep(stats),quote=FALSE)
  cat("\n")

  print(c(x$coef, x$icoef[2]), digits=digits)

  correl <- x$correl
  if (correlation && !is.null(x$correl))
    {  ## FEH
      p <- dim(correl)[2]
      if (p > 1)
        {
          cat("\nCorrelation of Coefficients:\n")
          ll <- lower.tri(correl)
          correl[ll] <- format(round(correl[ll], digits = digits))
          correl[!ll] <- ""
          print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
    }
  cat("\n")
  
  invisible()
}

vcov.pphsm <- function(object, ...) .NotYetImplemented()
