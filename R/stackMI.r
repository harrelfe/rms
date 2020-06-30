##' Bayesian Model Fitting and Stacking for Multiple Imputation
##'
##' Runs an \code{rms} package Bayesian fitting function such as \code{blrm} separately for each completed dataset given a multiple imputation result such as one produced by \code{Hmisc::aregImpute}.  Stacks the posterior draws and diagnostics across all imputations, and computes parameter summaries on the stacked posterior draws.
##' stackMI 
##' @param formula a model formula
##' @param fitter a Bayesian fitter
##' @param xtrans an object created by \code{transcan}, \code{\link{aregImpute}}, or \code{\link[mice]{mice}}
##' @param data data frame
##' @param n.impute number of imputations to run, default is the number saved in \code{xtrans}
##' @param dtrans see \code{Hmisc::fit.mult.impute}
##' @param derived see \code{Hmisc::fit.mult.impute}
##' @param subset an integer or logical vector specifying the subset of observations to fit
##' @param refresh see [rstan::sampling].  The default is 0, indicating that no progress notes are output.  If \code{refresh > 0} and \code{progress} is not \code{''}, progress output will be appended to file \code{progress}.  The default file name is \code{'stan-progress.txt'}.
##' @param progress see \code{refresh}.  Defaults to \code{''} if \code{refresh = 0}.  Note: If running interactively but not under RStudio, \code{rstan} will open a browser window for monitoring progress.
##' @param ... arguments passed to \code{fitter}
##' @return an \code{rms} fit object with expanded posterior draws and diagnostics
##' @author Frank Harrell
stackMI <-
  function(formula, fitter, xtrans, data,
           n.impute=xtrans$n.impute,
           dtrans, derived, subset,
           refresh=0,
           progress=if(refresh > 0) 'stan-progress.txt' else '', ...) {
    
  call <- match.call()

  using.Design <- FALSE
  used.mice <- any(class(xtrans)=='mids')
  if(used.mice && missing(n.impute)) n.impute <- xtrans$m

    draws       <- omega <- gammas <- eps <- NULL
    diagnostics <- list()
  
  for(i in 1 : n.impute) {
    if(progress != '')
      cat('\n-----------------------------------------------
Fitting imputed dataset number', i, 'of', n.impute, '\n\n',
                           file=progress, append=TRUE)

    if(used.mice) {
      completed.data <- mice::complete(xtrans, i)
      for(impvar in names(completed.data))
        if(length(attr(completed.data[[impvar]], 'contrasts')))
          attr(completed.data[[impvar]], 'contrasts') <- NULL
    }
    else {
      completed.data <- data
      imputed.data <-
        impute.transcan(xtrans, imputation=i, data=data,
                        list.out=TRUE, pr=FALSE, check=FALSE)
      ## impute.transcan works for aregImpute
      completed.data[names(imputed.data)] <- imputed.data
    }

    if(!missing(dtrans)) completed.data <- dtrans(completed.data)

    if(!missing(derived)) {
      stop('derived variables in fit.mult.imputed not yet implemented')
      eval(derived, completed.data)
    }

    if(using.Design) options(Design.attr=da)

    f <- if(missing(subset))
           fitter(formula, data=completed.data,
                  refresh=refresh, progress=progress, ...)
         else fitter(formula, data=completed.data[subset,],
                     refresh=refresh, progress=progress, ...)

    f$rstan <- NULL
    draws   <- rbind(draws, f$draws)
    omega   <- rbind(omega, f$omega)
    if(length(f$gammas))
      gammas <- if(i == 1) f$gammas else gammas + f$gammas
    if(length(f$eps))
      eps    <- if(i == 1) f$eps    else eps + f$eps

    dx <- f$diagnostics
    fail <- length(dx) && is.list(dx) && length(dx$fail) && dx$fail
    if(! fail) diagnostics[[i]] <- dx

    if(inherits(f, 'rms')) {
      using.Design <- TRUE
      da           <- f$Design
      }
   }

    if(length(gammas)) gammas <- gammas / n.impute
    if(length(eps))    eps    <- eps    / n.impute
    
    param <- rbind(mean=colMeans(draws), median=apply(draws, 2, median))

    f$n.impute    <- n.impute
    f$param       <- param
    f$draws       <- draws
    f$omega       <- omega
    f$gammas      <- gammas
    f$eps         <- eps
    f$diagnostics <- diagnostics

    f$formula <- formula
#    f$assign  <- assign
    f$call    <- call
    if(using.Design) options(Design.attr=NULL)
    class(f) <- c('stackImpute', class(f))
    f
  }
