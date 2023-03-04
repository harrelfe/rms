##' Process Special Multiple Imputation Output
##'
##' Processes lists that have one element per imputation
##' @title processMI
##' @param object a fit object created by [Hmisc::fit.mult.impute()]
##' @param ... ignored
##' @return an object that resembles something created by a single fit without multiple imputation
##' @seealso [processMI.fit.mult.impute()]
##' @author Frank Harrell
##' @md
processMI <- function(object, ...) UseMethod("processMI")

##' Process Special Multiple Imputation Output From `fit.mult.impute`
##'
##' Processes a `funresults` object stored in a fit object created by `fit.mult.impute` when its `fun` argument was used.  These objects are typically named `validate` or `calibrate` and represent bootstrap or cross-validations run separately for each imputation.  See [this](https://hbiostat.org/rmsc/validate.html#sec-val-mival) for a case study.
##' @title processMI.fit.mult.impute
##' @param object a fit object created by `fit.mult.impute`
##' @param which specifies which component of the extra output should be processed
##' @param plotall set to `FALSE` when `which='calibrate'` to suppress having `ggplot` render a graph showing calibration curves produced separately for all the imputations
##' @param nind set to a positive integer to use base graphics to plot a matrix of graphs, one each for the first `nind` imputations, and the overall average calibration curve at the end
##' @param ... ignored
##' @return an object like a `validate` or `calibrate` result obtained when no multiple imputation was done.  This object is suitable for `print` and `plot` methods for these kinds of resampling validation objects.
##' @seealso [Hmisc::fit.mult.impute()]
##' @md
##' @author Frank Harrell
processMI.fit.mult.impute <-
  function(object, which=c('validate', 'calibrate'),
           plotall=TRUE, nind=0, ...) {
    which <- match.arg(which)

    r <- lapply(object$funresults,
                function(x)
                  if(which %nin% names(x))
                    stop(paste('fun result does not contain', which))
                  else x[[which]] )
    
    n.impute <- length(r)
    
    if(which == 'validate') {
      v <- r[[1]]
      if(n.impute > 1) for(i in 2 : n.impute) v <- v + r[[i]]
      ## Average all indexes but sum the number of resamples
      v[, colnames(v) != 'n'] <- v[, colnames(v) != 'n'] / n.impute
      attr(v, 'n.impute') <- n.impute
      return(v)
      }
    if(which == 'calibrate') {
      cal  <- inherits(r[[1]], 'calibrate')
      cald <- inherits(r[[1]], 'calibrate.default')
      km   <- cal && ('KM' %in% colnames(r[[1]]))
      predname <- if(cald) 'predy'
                  else if(cal && km) 'mean.predicted'
                  else 'pred'
      if(cald || (cal && TRUE)) {
        ## Create a tall data frame with all predicted values and overfitting-
        ## corrected estimates
        d <- NULL
        for(i in 1 : n.impute)
          d <- rbind(d,
                     data.frame(imputation=i,
                                predicted=r[[i]][, predname],
                                corrected=r[[i]][, if(km) 'KM.corrected'
                                                   else 'calibrated.corrected']) )
        if(plotall) {
          g <- ggplot(d, aes(x=predicted, y=corrected,
                             color=factor(imputation))) +
            geom_line(alpha=0.3) + xlab('Predicted') +
            ylab('Estimated Actual, Overfitting-Corrected') +
            guides(color=FALSE)
          print(g)
          }
        ## Find range of predicted values over all imputations
        np <- nrow(r[[1]])
        ## Compute new common grid to interpolate all imputations to
        if(km) {
          pred <- sort(unique(d$predicted))
          ## Remove all points that are within than 0.005 to previous
          ## by rounding to nearest 0.005
          pred <- unique(round(pred / 0.005) * 0.005)
          }
        else pred <- seq(min(d$predicted), max(d$predicted), length=np)
        np <- length(pred)
        ## For each imputation interpolate all values to this grid
        ## Accumulate sums of interpolated values and get the final
        ## result by averaging
        k <- matrix(0, nrow=np, ncol=ncol(r[[1]]),
                    dimnames=list(NULL, colnames(r[[1]])))
        k[, predname] <- pred
        B <- 0
        for(i in 1 : n.impute) {
          x <- r[[i]]
          B <- B + attr(x, 'B')
          for(j in setdiff(colnames(k), predname))
            k[, j] <- k[, j] + approxExtrap(x[, predname], x[, j],
                                            xout=pred)$y
        }
        for(j in setdiff(colnames(k), c(predname, 'n')))
          k[, j] <- k[, j] / n.impute
        at <- attributes(r[[1]])
        at$dim <- at$dimnames <- NULL
        attributes(k) <- c(attributes(k), at)
        attr(k, 'B')  <- B
        if(nind > 0) {
          oldpar <- par(mfrow=grDevices::n2mfrow(nind + 1))
          on.exit(par(oldpar))
          for(i in 1 : nind)
            plot(r[[i]], main=paste('Imputation', i))
            plot(k,      main=paste('Average Over', n.impute, 'Imputations'))
          }
        return(k)
      }
      stop(paste('calibrate object class', paste(class(r[[1]]), collapse=' '),
                 'not yet implemented'))
    }
  }

utils::globalVariables(c('predicted', 'corrected', 'imputation'))
