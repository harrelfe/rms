##' Relative Explained Variation
##'
##' Computes measures of relative explained variation for each predictor in an `rms` or `rmsb` model fit `object`.  This is similar to `plot(anova(fit), what='proportion R2')`.  For an `ols` model the result is exactly that.  Uncertainty intervals are computed if the model fit is from `rmsb` or was run through [bootcov()] with `coef.reps=TRUE`.  The results may be printed, and there is also a `plot` method.
##'
##' When `object` is not an `ols` fit, the linear predictor from the fit in `object` is predicted from the original predictors, resulting in a linear model with \eqn{R^{2}=1.0}.  The partial \eqn{R^2} for each predictor from a new `ols` fit is the relative explained variation.  The process is repeated when bootstrap coefficients repetitions or posterior draws are present, to get uncertainty intervals.  So relative explained variation is the proportion of variation in the initial model's predicted values (on the linear predictor scale) that is due to each predictor.
##'
##' Nonlinear and interaction terms are pooled with main linear effect of predictors, so relative explained variation for a predictor measures its total impact on predicted values, either as main effects or effect modifiers (interaction components).
##' @title rexVar
##' @param object a fit from `rms` or `rmsb` 
##' @param data a data frame, data table, or list providing the predictors used in the original fit
##' @param ns maximum number of bootstrap repetitions or posterior draws to use
##' @param cint confidence interval coverage probability for nonparametric bootstrap percentile intervals, or probability for a Bayesian highest posterior density interval for the relative explained variations.
##' @return a vector (if bootstrapping or Bayesian posterior sampling was not done) or a matrix otherwise, with rows corresponding to predictors and colums `REV`, `Lower`, `Upper`.  The returned object is of class `rexVar`.
##' @author Frank Harrell
##' @md
##' @examples
##' set.seed(1)
##' n <- 100
##' x1 <- rnorm(n)
##' x2 <- rnorm(n)
##' x3 <- rnorm(n)
##' y  <- x1 + x2 + rnorm(n) / 2.
##' d <- data.frame(x1, x2, x3, y)
##' dd <- datadist(d); options(datadist='dd')
##' f  <- ols(y ~ pol(x1, 2) * pol(x2, 2) + x3,
##'           data=d, x=TRUE, y=TRUE)
##' plot(anova(f), what='proportion R2', pl=FALSE)
##' rexVar(f)
##' g <- bootcov(f, B=20, coef.reps=TRUE)
##' rexVar(g, data=d)
##' f <- orm(y ~ pol(x1,2) * pol(x2, 2) + x3,
##'          data=d, x=TRUE, y=TRUE)
##' rexVar(f, data=d)
##' g <- bootcov(f, B=20, coef.reps=TRUE)
##' rexVar(g, data=d)
##' \dontrun{
##' require(rmsb)
##' h <- blrm(y ~ pol(x1,2) * pol(x2, 2) + x3, data=d)
##' rexVar(h, data=d)
##' }
##' options(datadist=NULL)

rexVar <- function(object, data, ns=500, cint=0.95) {
  rex <- function(olsfit) {
    a <- anova(olsfit)
    
    pss <- a[, 'Partial SS']
    ssr <- a['TOTAL', 'Partial SS']
    sst <- a['ERROR', 'Partial SS'] + ssr

    rm <- c("TOTAL NONLINEAR","TOTAL NONLINEAR + INTERACTION",
            "TOTAL INTERACTION","TOTAL", 
            " Nonlinear"," All Interactions", "ERROR",
            " f(A,B) vs. Af(B) + Bg(A)")
    rn <- rownames(a)
    rm <- c(rm, rn[substring(rn, 2, 10) == "Nonlinear"])
    pss <- a[rn %nin% rm, 'Partial SS']
    names(pss) <- sub(' (Factor+Higher Order Factors)', '', names(pss),
                      fixed=TRUE)
    r <- pss / ssr
    names(r) <- trimws(names(r))
    r
    }
    
  draws <- object$draws
  drawtype <- 'bayes'
  if(! length(draws)) {
    draws <- object$boot.Coef
    drawtype <- 'bootstrap'
    }
  if(inherits(object, 'ols') && ! length(draws))
    return(structure(rex(object), class='rexVar'))

  .lp.  <- if(inherits(object, 'blrm')) predict(object, cint=FALSE)
  else predict(object)
  form  <- formula(object)
  form  <- update(form, .lp. ~ .)   # replace dependent variable with .lp.
  data$.lp. <- .lp.
  f     <- ols(form, data=data, x=TRUE)
  X     <- f$x
  overall.rex <- rex(f)
  ## If not a Bayesian fit from rmsb or bootstrapped by robcov, finished
  if(! length(draws)) return(structure(overall.rex, class='rexVar'))

  ni <- num.intercepts(object)

  ## Recreate linear predictor for posterior draws or bootstraps
  ns <- min(ns, nrow(draws))
  rx <- matrix(NA, nrow=ns, ncol=length(overall.rex))
  for(i in 1 : ns) {
    ## matxv drops excess coefficients representing intercepts not
    ## used in X
    .lp. <- matxv(X, draws[i, , drop=TRUE])
    g <- lm.fit.qr.bare(X, as.vector(.lp.),
                        tolerance=1e-13, intercept=TRUE, xpxi=TRUE)
    ## Trick to update original ols model fit; speeds up over
    ## having to recreate design matrix each time
    f$coefficients   <- g$coefficients
    f$var            <- g$xpxi
    f$stats['Sigma'] <- 1.0   # really zero
    rx[i, ] <- rex(f)
  }

  lim <- apply(rx, 2, rmsb::HPDint, prob=cint)
  r   <- cbind(REV=overall.rex, Lower=lim[1, ], Upper=lim[2, ])
  rownames(r) <- names(overall.rex)
  structure(r, class='rexVar', drawtype=drawtype)
}

##' Print rexVar Result
##'
##' Prints the results of an `rexVar` call
##' @title print.rexVar
##' @param x a vector or matrix created by `rexVar`
##' @param title character string which can be set to `NULL` or `''` to suppress
##' @param digits passed to [round()]
##' @param ... unused
##' @return invisible
##' @author Frank Harrell
##' @md
print.rexVar <- function(x, title='Relative Explained Variation',
                         digits=3, ...) {
  if(length(title) && title != '')
    cat('\nRelative Explained Variation\n\n')
  attr(x, 'drawtype') <- NULL
  invisible(print(round(unclass(x), digits)))
  }

##' Plot rexVar Result
##'
##' Makes a dot chart displaying the results of `rexVar`.  Base graphics are used unless `options(grType='plotly')` is in effect, in which case a `plotly` graphic is produced with hovertext
##' @title plot.rexVar
##' @param x a vector or matrix created by `rexVar`
##' @param xlab x-axis label
##' @param xlim x-axis limits; defaults to range of all values (limits and point estimates)
##' @param pch plotting symbol for dot
##' @param sort defaults to sorted predictors in descending order of relative explained variable.  Can set to `ascending` or `none`.
##' @param margin set to `TRUE` to show the REV values in the right margin if using base graphics
##' @param height optional height in pixels for `plotly` graph 
##' @param width likewise optional width
##' @param ... arguments passed to `dotchart2` or `dotchartpl`
##' @return `plotly` graphics object if using `plotly`
##' @author Frank Harrell
##' @md
plot.rexVar <- function(x, xlab='Relative Explained Variation',
                        xlim=NULL,
                        pch=16, sort=c("descending", "ascending", "none"),
                        margin=FALSE,  height=NULL, width=NULL, ...) {
    
  sort <- match.arg(sort)
  isbase <- Hmisc::grType() == 'base'
  if(! is.matrix(x))
    x <- matrix(x, ncol=1, dimnames=list(names(x), 'REV'))
  nr <- nrow(x)
  if(! isbase && ! length(height))
    height <- plotlyParm$heightDotchart(nr)
  drawtype <- attr(x, 'drawtype')

  i <- switch(sort,
              none       = 1 : nr,
              descending = order(x[, 'REV'], decreasing=TRUE),
              ascending  = order(x[, 'REV']))
  x <- x[i,, drop=FALSE]
  rownames(x) <- trimws(rownames(x))
  if(! length(xlim)) xlim <- range(x)
  ul <- ncol(x) > 1

  if(isbase) {
    if(margin)
    dotchart2(as.vector(x[, 'REV']), labels=rownames(x),
              xlab=xlab, pch=pch, xlim=xlim,
              auxdata=if(margin) round(x[, 'REV'], 3),
              ...)
    else dotchart2(as.vector(x[, 'REV']), labels=rownames(x),
                   xlab=xlab, pch=pch, xlim=xlim,
                   ...)

    if(ul) {
      dotchart2(x[, 'Lower'], pch=91, add=TRUE)
      dotchart2(x[, 'Upper'], pch=93, add=TRUE)
      }
    return(invisible())
  }
  
  lname <- if(length(drawtype))
             switch(drawtype,
                    bayes='HPD Interval',
                    bootstrap='Bootstrap CI')
  
  dotchartpl(x[, 'REV'], major=rownames(x),
             lower=if(ul) x[,'Lower'],
             upper=if(ul) x[,'Upper'],
             htext=format(round(x[, 'REV'], 3)),
             xlab=xlab, xlim=xlim,
             limitstracename=lname,
             width=width, height=height, ...)
}
