##' Update Model LR Statistics After Multiple Imputation
##'
##' For fits from `orm, lrm, orm, cph, psm` that were created using `fit.mult.impute` with `lrt=TRUE` or equivalent options and for which `anova` was obtained using `processMI(fit, 'anova')` to compute imputation-adjusted LR statistics.  `LRupdate` uses the last line of the `anova` result (containing the overall model LR chi-square) to update `Model L.R.` in the fit `stats` component, and to adjust any of the new R-square measures in `stats`.
##' 
##' For models using Nagelkerke's R-squared, these are set to `NA` as they would need to be recomputed with a new intercept-only log-likelihood, which is not computed by `anova`.  For `ols` models, R-squared is left alone as it is sample-size-independent and `print.ols` prints the correct adjusted R-squared due to `fit.mult.impute` correcting the residual d.f. in stacked fits.
##' @title LRupdate
##' @param fit an `rms` fit object
##' @param anova the result of `processMI(..., 'anova')`
##' @return new fit object like `fit` but with the substitutions made
##' @author Frank Harrell
##' @seealso [processMI.fit.mult.impute()], [Hmisc::R2Measures()]
##' @md
##' @examples
##' \dontrun{
##' a <- aregImpute(~ y + x1 + x2, n.impute=30, data=d)
##' f <- fit.mult.impute(y ~ x1 + x2, lrm, a, data=d, lrt=TRUE)
##' a <- processMI(f, 'anova')
##' f <- LRupdate(f, a)
##' print(f, r2=1:4)   # print all imputation-corrected R2 measures
##' }
LRupdate <- function(fit, anova) {
  cl <- class(fit)
  if('rms' %nin% cl) stop('fit is not an rms fit')
  if(! inherits(anova, 'anova.rms'))
    stop('anova is not the result of rms anova()')

  lr <- anova['TOTAL', 'Chi-Square']
  s  <- fit$stats
  s['Model L.R.'] <- lr
  df <- s['d.f.']
  s['P'] <- pchisq(lr, df, lower.tail=FALSE)
  n  <- if('Obs' %in% names(s)) s['Obs'] else s['n']  # n for ols
  mod <- ''
  for(j in c('ols', 'lrm', 'orm', 'cph', 'psm')) if(j %in% cl) mod <- j
  if(mod == '')
    stop('model not from ols, lrm, orm, cph, psm')

  if(mod != 'ols') {
    if(mod %in% c('lrm', 'orm', 'cph', 'psm') && 'R2' %in% names(s))
      s['R2'] <- NA    # Nagelkerke R^2 no longer correct

    ## Fetch effective sample size (scalar) or Y frequency distribution
    ess <- switch(mod,
                  lrm = fit$freq,
                  orm = fit$freq,
                  cph = s['Events'],
                  psm = s['Events'])

    r2m <- R2Measures(lr, df, n, ess)
    i <- grep('R2\\(', names(s))
    if(length(i) != length(r2m))
      stop('number of R2 from R2Measures (', length(r2m),
           ') does not equal number stored in fit (', length(i), ')')
    s[i] <- r2m
    names(s)[i] <- names(r2m)
  }
  
  fit$stats <- s
  fit
  }
