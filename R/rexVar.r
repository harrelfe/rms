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
    names(pss) <- sub(' (Factor+Higher Order Factors)', '', names(pss), fixed=TRUE)
    pss / ssr
    }
    
  draws <- object$draws
  if(! length(draws)) draws <- object$boot.Coef
  if(inherits(object, 'ols') && ! length(draws)) return(rex(object))

  .lp.  <- predict(object)
  if(is.list(.lp.)) .lp. <- .lp.$linear.predictors  # for blrm fits
  form  <- formula(object)
  form  <- update(form, .lp. ~ .)
  data$.lp. <- .lp.
  f     <- ols(form, data=data, x=TRUE)
  X     <- f$x
  overall.rex <- rex(f)
  ## If not a Bayesian fit from rmsb or bootstrapped by robcov, finished
  if(! length(draws)) return(overall.rex)

  if((ncol(draws) != length(coef(object))) ||
     ! identical(colnames(draws), names(coef(object))))
    stop('coefficients in fit object are not compatible with columns of posterior draws or bootstraps')

  ## Recreate linear predictor for posterior draws or bootstraps
  ns <- min(ns, nrow(draws))
  rx <- matrix(NA, nrow=ns, ncol=length(overall.rex))
  for(i in 1 : ns) {
    .lp. <- matxv(X, draws[i, ])
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
  r
}
