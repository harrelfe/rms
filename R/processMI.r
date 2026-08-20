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
##'
##' For `which='anova'` assumes that the `fun` given to `fit.mult.impute` runs `anova(fit, test='LR')` to get likelihood ratio tests, and that `method='stack'` was specified also so that a final `anova` was run on the stacked combination of all completed datasets.  The method of [Chan and Meng (2022)](https://hbiostat.org/rmsc/missing.html#sec-missing-lrt) is used to obtain overall likelihood ratio tests, with each line of the `anova` table getting a customized adjustment based on the amount of missing information pertaining to the variables tested in that line.  The resulting statistics are chi-square and not $F$ statistics as used by Chan and Meng.  This will matter when the estimated denominator degrees of freedom for a variable is small (e.g., less than 50).  These d.f. are reported so that user can take appropriate cautions such as increasing `n.impute` for `aregImpute`.
##' @title processMI.fit.mult.impute
##' @param object a fit object created by `fit.mult.impute`
##' @param which specifies which component of the extra output should be processed
##' @param plotall set to `FALSE` when `which='calibrate'` to suppress having `ggplot` render a graph showing calibration curves produced separately for all the imputations
##' @param nind set to a positive integer to use base graphics to plot a matrix of graphs, one each for the first `nind` imputations, and the overall average calibration curve at the end
##' @param prmi set to `FALSE` to not print imputation corrections for `anova`
##' @param ... ignored
##' @return an object like a `validate`, `calibrate`, or `anova` result obtained when no multiple imputation was done.  This object is suitable for `print` and `plot` methods for these kinds of objects.
##' @seealso [Hmisc::fit.mult.impute()]
##' @md
##' @author Frank Harrell
processMI.fit.mult.impute <-
  function(object, which=c('validate', 'calibrate', 'anova'),
           plotall=TRUE, nind=0, prmi=TRUE, ...) {
    which <- match.arg(which)

    r <- lapply(object$funresults,
                function(x)
                  if(which %nin% names(x))
                    stop(paste('fun result does not contain', which))
                  else x[[which]] )
    
    n.impute <- object$n.impute
    if(which == 'anova' && length(r) != (n.impute + 1))
      stop('runresults has wrong length for anova')
    if(which != 'anova' && length(r) != n.impute)
      stop('runresults has wrong length for non-anova')
    
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
        else pred <- seq(min(d$predicted), max(d$predicted), length.out=np)
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
    if(which == 'anova') {
      M <- n.impute
      ## Get number of tests done by anova
      nt <- nrow(r[[1]])
      ## Compute mean of all n.impute LR chi-squares
      lrt <- rep(0., nt)
      for(i in 1 : M) lrt <- lrt + r[[i]][, 'Chi-Square']
      lrt <- lrt / M
      ## Get result from stack dataset model fit
      A   <- r[[n.impute + 1]]
      LRT <- A[, 'Chi-Square'] / M
      df  <- r[[n.impute + 1]][, 'd.f.']
      ## For each test do the MI corrections
      rhat <- pmax(0., ((M + 1) / (df * (M - 1))) * (lrt - LRT))
      fhat <- rhat / (1. + rhat)
      df2  <- df * (M - 1) / fhat / fhat
      dfact <- 1. / (1. + rhat)
      mi.info <- data.frame(Test = rownames(A),
                            'Missing Information' = fhat,
                            'Denominator d.f.'    = df2,
                            'Chi-Square Discount' = dfact,
                            check.names=FALSE)
      attr(A, 'mi.info') <- mi.info
      A[, 'Chi-Square'] <- dfact * LRT
      A[, 'P']          <- pchisq(dfact * LRT, df, lower.tail=FALSE)

      A
    }
  }

##' Print Information About Impact of Imputation
##'
##' For the results of `processMI.fit.mult.impute` prints or renders a formatted table (html, LaTeX, or Typst, according to `options(prType=)`) summarizing various correction factors related to missing data multiple imputation.
##' @title prmiInfo
##' @param x an object created by `processMI(..., 'anova')`
##' @return nothing 
##' @author Frank Harrell
##' @md
##' @examples
##' \dontrun{
##' a <- aregImpute(...)
##' f <- fit.mult.impute(...)
##' v <- processMI(f, 'anova')
##' prmiInfo(v)
##' }
prmiInfo <- function(x) {
  m <- attr(x, 'mi.info')
  if(! length(m)) stop('object does not have mi.info attributes')
  
  for(j in 2:4) m[, j] <- format(round(m[, j], c(NA,3,1,3)[j]))

  lang <- prType()
  if(lang == 'plain') { cat('\n'); print(m); cat('\n'); return(invisible()) }
  ## unified latex/html/typst rendering -- see end-of-file note

  specs <- markupSpecs[[lang]]
  rowl  <- m$Test
  if('MS' %in% names(m)) rowl[rowl=='TOTAL'] <- 'REGRESSION'
  bold  <- specs$bold

  ## row-label translation, extended to latex/typst
  rowl <- switch(lang, latex = latexTranslate(rowl),
                 typst = typstTranslate(rowl), rowl)

  ## Translate interaction symbol (*) to times symbol. Searches for the
  ## already-escaped pattern under typst, matching anova.rms's identical
  ## fix (typstTranslate() above has already escaped bare * to \*).
  pat  <- if(lang == 'typst') '\\*' else '*'
  rowl <- gsub(pat, specs$times, rowl, fixed=TRUE)

  ## Put TOTAL rows in boldface
  rowl <- ifelse(substring(rowl, 1, 5) %in% c("REGRE", "ERROR", "TOTAL"),
                 bold(rowl), rowl)
  rowl <- ifelse(substring(rowl, 1, 1) == " ",
               paste0(specs$lspace, specs$italics(substring(rowl,2)), sep=""),
               rowl) # preserve leading blank
  m$Test <- rowl

  ## chisq() wrapping is format-inconsistent (typst pre-wrapped, latex
  ## bare) -- same fix as anova.rms's chi-square column name; html's
  ## chisq() is already display-ready and untouched, matching original
  chi2lab <- specs$chisq(add='')
  if(lang != 'html' && substring(chi2lab, 1, 1) != '$')
    chi2lab <- specs$math(chi2lab)

  names(m) <- c('Test', 'Missing\nInformation\nFraction',
                'Denominator\nd.f.', paste(chi2lab, 'Discount'))

  rownames(m) <- m$Test
  m$Test <- NULL

  w <- rms_tiny_table(m, lang, caption = 'Imputation penalties')
  if(lang == 'html') rendHTML(w) else rendHTML(w, html = FALSE)
  }


utils::globalVariables(c('predicted', 'corrected', 'imputation'))

## -----------------------------------------------------------------------
## Detailed notes on extending prmiInfo's table output to latex/typst,
## and unifying it with rms_tiny_table, referenced by the short comment
## inline above. Most of this is simply building out latex/typst
## coverage for the first time -- part of this project's broader goal,
## not corrections to anything broken in the previously html-only
## version.
##
## - prmiInfo was previously html-only (if(prType()=='html') {...} else
##   {plain-text print}) -- no latex or typst branch at all, meaning
##   fit.mult.impute imputation-penalty tables never rendered as a
##   formatted table under those two prTypes, falling straight through
##   to plain-text printing regardless. Restructured to the same
##   simplified dispatch pattern already used by print.anova.rms/
##   print.summary.rms/print.validate after their own refactor
##   (if(lang == 'plain') {...} else {...unified rendering...}), rather
##   than a three-way html/latex/typst switch.
##
## - Row-label handling (interaction-symbol substitution, bold TOTAL
##   rows, italicized/indented sub-rows) is structurally identical to
##   anova.rms's own latex.anova.rms() -- this function appears to have
##   been directly modeled on it. Extended using the exact same fixes
##   already established and compile-tested there, rather than
##   reinventing them: latexTranslate()/typstTranslate() applied to
##   rowl before the interaction-symbol substitution (previously absent
##   entirely, since only html was ever supported and html's specs
##   didn't need it); the substitution itself searches for the
##   already-escaped "\\*" pattern under typst rather than a bare "*"
##   (typstTranslate() has already escaped bare * by that point, since *
##   is Typst's bold/emphasis shorthand -- searching for the bare
##   pattern would match the * inside \*, leaving an orphaned backslash,
##   the same bug class found and fixed in anova.rms and reused
##   verbatim here).
##
## - chi-square column label: specs$chisq()'s $-wrapping is
##   format-inconsistent (typst pre-wrapped, latex bare, html
##   already display-ready and not math-mode at all) -- same issue
##   and same fix already established in anova.rms's note [5]: check
##   whether the returned value already starts with "$" rather than
##   assume either format's convention, wrapping with specs$math()
##   only when it doesn't. html is left completely untouched, matching
##   its original behavior exactly (specs$chisq() used directly, no
##   wrapping, since it was already correct as html markup).
##
## - Rendering: htmlTable::htmlTable() (html-only, hand-built
##   css.table/css.cell font-size and padding strings) replaced with
##   rms_tiny_table(m, lang, caption='Imputation penalties') for all
##   three formats -- no new rms_tiny_table parameters were needed here;
##   this table's shape (row labels + right-aligned numeric columns +
##   caption) is exactly what rms_tiny_table already supported by
##   default; the row-label/data-column alignment split ('l' for Test,
##   'r' for the three statistic columns) matches rms_tiny_table's
##   existing defaults exactly, so no align_data override was needed
##   either.
##
## - Multi-line column headers ("Missing\nInformation\nFraction", etc.)
##   use "\n" uniformly across all three formats, including html --
##   matching summary.rms's already-established, working pattern for
##   this, rather than the original code's own html-specific "<br>"
##   convention. rms_tiny_table_core's unconditional
##   format_tt(i=0, linebreak="\n") handles "\n" correctly for html too,
##   making a separate <br>-based header string unnecessary.
## -----------------------------------------------------------------------
