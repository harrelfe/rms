validate <-
  function(fit,  method="boot", B=40,
           bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0,
           force=NULL, estimates=TRUE, pr=FALSE,...)
  UseMethod("validate")


#Resampling optimism of discrimination and reliability of an ols regression
#B: # reps
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each bootstrap rep
#Requires: predab.resample, fastbw, ols
#Frank Harrell 11 June 91

validate.ols <- function(fit, method="boot",
	 B=40, bw=FALSE, rule="aic", type="residual",
	 sls=.05, aics=0, force=NULL, estimates=TRUE,
    pr=FALSE, u=NULL, rel=">", tolerance=1e-7, ...)
{
  fit.orig <- fit

  penalty.matrix <- fit.orig$penalty.matrix

  discrim <- function(x, y, fit, iter, evalfit=FALSE, u=NULL, rel=NULL,
                      pr=FALSE, ...)
	{
      resid <- if(evalfit) fit$residuals else y - x
      n <- length(resid)
      sst <- (n - 1) * var(y)   # sum(y^2) - (sum(y)^2)/n
      mse <- sum(resid^2)
      rsquare <- 1 - mse / sst
      mse <- mse / n

      if(evalfit) {	#Fit being examined on sample used to fit
        intercept <- 0
        slope     <- 1
      } else {
        if(length(fit$coef)==1) {intercept <- mean(y)-mean(x); slope <- 1}
        else {
          coef <- lsfit(x, y)$coef   #Note x is really x*beta from other fit
          intercept <- coef[1]
          slope     <- coef[2]
        }
      }

      z <- c(rsquare, mse, GiniMd(slope*x), intercept, slope)
      nam <- c("R-square", "MSE", "g", "Intercept", "Slope")
      if(length(u)) {
        yy <- if(rel==">") ifelse(y >  u, 1, 0)
        else if(rel==">=") ifelse(y >= u, 1, 0)
        else if(rel=="<")  ifelse(y <  u, 1, 0)
        else ifelse(y <= u, 1, 0)
        z <- c(z, somers2(x,yy)["Dxy"])
        nam <- c(nam, paste("Dxy Y", rel, format(u), sep=""))
        if(rel == ">" | rel == ">=") P <- pnorm(- (u - x) / sqrt(mse))
        else P <- pnorm((u - x) / sqrt(mse))
        P0  <- sum(yy) / n
        L   <- -2 * sum(yy * logb(P)  + (1 - yy) * logb(1 - P ))
        L0  <- -2 *sum(yy * logb(P0) + (1 - yy) * logb(1 - P0))
        R2  <- (1 - exp(-(L0 - L) / n)) / (1 - exp(-L0 / n))
        z   <- c(z, R2)
        nam <- c(nam, paste("R2 Y", rel, format(u), sep=""))
      }
      names(z) <- nam
      z
    }

  ols.fit <- function(x, y, tolerance=1e-7, backward,
                      penalty.matrix=NULL, xcol=NULL, ...)
  {
    fail <- FALSE
    if(!length(x)) {
      ybar <- mean(y)
      n <- length(y)
      residuals <- y - ybar
      v <- sum(residuals ^ 2) / (n - 1)
      return(list(coef=ybar, var=v / n, residuals=residuals, fail=fail))
    }
    if(length(penalty.matrix) > 0) {
      if(length(xcol)) {
        xcol <- xcol[-1] - 1   # remove position for intercept
        penalty.matrix <- penalty.matrix[xcol, xcol, drop=FALSE]
      }
      fit <- lm.pfit(x, y, penalty.matrix=penalty.matrix,
                     tol=tolerance)
      if(any(is.na(fit$coefficients))) fail <- TRUE
    }
    else {
      fit <- lm.fit.qr.bare(x, as.vector(y), tolerance=tolerance,
                            intercept=TRUE, xpxi=TRUE)
      if(any(is.na(fit$coefficients))) fail <- TRUE
      if(backward)
        fit$var <- sum(fit$residuals^2) * fit$xpxi/
          (length(y) - length(fit$coefficients))
    }
    c(fit, fail=fail)
  }

  predab.resample(fit.orig, method=method, fit=ols.fit, measure=discrim,
                  pr=pr, B=B, bw=bw, rule=rule, type=type, sls=sls, aics=aics,
                  force=force, estimates=estimates, tolerance=tolerance,
                  backward=bw,u=u, penalty.matrix=penalty.matrix,
                  rel=rel, ...)
}

print.validate <- function(x, digits=4, B=Inf, ...)
{

  lang <- prType()
  ## unified dispatch for latex/html/typst -- see note [1] at end of
  ## file for why this collapsed from a two-way check
  if(lang != 'plain') return(latex.validate(x, digits=digits, B=B, ...))

    kept <- attr(x, 'kept'); attr(x, 'kept') <- NULL
    print(round(unclass(x), digits), ...)
    if(length(kept) && B > 0) {
      cat("\nFactors Retained in Backwards Elimination\n\n")
      varin <- ifelse(kept, '*', ' ')
      print(varin[1:min(nrow(varin), B),], quote=FALSE)
      cat("\nFrequencies of Numbers of Factors Retained\n\n")
      nkept <- apply(kept, 1, sum)
      tkept <- table(nkept)
      names(dimnames(tkept)) <- NULL
      print(tkept)
    }
  }

## latex.validate: builds and renders up to three tables (main
## validation statistics; factors-retained marker table; frequencies of
## factors retained) for the current output language (latex/html/
## typst). Basic calling sequence unchanged from the original
## latex.validate/html.validate -- html.validate below now delegates
## here, matching the pattern already used by html.anova.rms and
## html.summary.rms. See notes [2]-[5] at end of file for the reasoning
## behind unifying what were previously two separate, duplicated
## functions, the row-name translation additions, and the vestigial
## parameters kept for signature compatibility.
latex.validate <- function(object, digits=4, B=Inf, file='', append=FALSE,
                           title=first.word(deparse(substitute(x))),
                           caption=NULL, table.env=FALSE,
                           size='normalsize',
                           extracolsize=size, ...)
  {
    lang  <- prType()

    chg <- function(x, old, new) {
      names(new) <- old
      tx <- new[x]
      ifelse(is.na(tx), x, tx)
    }
    x <- object
    kept <- attr(x, 'kept'); attr(x, 'kept') <- NULL
    cn <- colnames(x)
    cn <- chg(cn, c('index.orig', 'training', 'test', 'optimism',
                    'index.corrected', 'Lower', 'Upper', 'n'),
              c('Original\nSample', 'Training\nSample',
                'Test\nSample', 'Optimism', 'Corrected\nIndex', 'Lower', 'Upper',
                'Successful\nResamples'))
    rn <- rownames(x)
    ## row-name translation, format-specific -- see note [3] at end of file
    rnNew <- switch(lang,
                    latex = c('$D_{xy}$','$R^{2}$','$R^{2}$', '$E_{\\max}$',
                              '$D$','$U$','$Q$','$B$','$g$','$g_{p}$',
                              '$g_{r}$','$\\rho$',
                              '$|\\overline{P(Y\\geq Y_{0.5})-\\frac{1}{2}}|$'),
                    html  = c('<i>D<sub>xy</sub></i>',
                              '<i>R<sup>2</sup></i>','<i>R<sup>2</sup></i>',
                              '<i>E<sub>max</sub></i>','<i>D</i>','<i>U</i>',
                              '<i>Q</i>','<i>B</i>','<i>g</i>',
                              '<i>g<sub>p</sub></i>',
                              '<i>g<sub>r</sub></i>','<i>&rho;</i>',
                              'Mean |P(Y&ge;Y<sub>0.5</sub>)-0.5|'),
                    typst = c('$D_("xy")$', '$R^(2)$', '$R^(2)$',
                              '$E_("max")$', '$D$', '$U$', '$Q$', '$B$',
                              '$g$', '$g_(p)$', '$g_(r)$', '$rho$',
                              '$|overline("P"(Y >= Y_(0.5)) - 1/2)|$'))
    rn <- chg(rn, c('Dxy','R2','R-square','Emax','D','U','Q','B','g','gp','gr',
                    'rho','pdm'), rnNew)
    dimnames(x) <- list(rn, cn)

    cdec <- ifelse(cn == 'Successful\nResamples', 0, digits)
    z <- unclass(x)
    ## Bug in htmlTable::txtRound for vector digits -- rounded directly
    ## here instead, format-agnostic
    for(i in 1 : length(cdec)) z[, i] <- round(z[, i], cdec[i])
    z <- as.data.frame(z, check.names = FALSE)
    rownames(z) <- rn

    R <- rms_tiny_table(z, lang, caption = caption)

    if(length(kept) && B > 0) {
      ## Bullet marker for retained variables, format-specific -- see
      ## note [4] at end of file
      bullet <- switch(lang, latex = '$\\bullet$',
                       html  = htmlSpecial('mediumsmallwhitecircle'),
                       typst = '\u2022')
      varin <- ifelse(kept, bullet, ' ')
      nr <- nrow(varin)
      varin <- varin[1:min(nrow(varin), B),, drop=FALSE]
      cap <- 'Factors Retained in Backwards Elimination'
      if(nr > B)
        cap <- switch(lang,
                      html = paste0(cap, '<br>First ', B, ' Resamples'),
                      paste0(cap, ' (first ', B, ' resamples)'))
      varin <- as.data.frame(varin, check.names = FALSE)
      R <- c(R, rms_tiny_table(varin, lang, caption = cap,
                               show_rownames = FALSE))

      cap <- 'Frequencies of Numbers of Factors Retained'
      nkept <- apply(kept, 1, sum)
      tkept <- t(as.matrix(table(nkept)))
      tkept <- as.data.frame(tkept, check.names = FALSE)
      R <- c(R, rms_tiny_table(tkept, lang, caption = cap,
                               show_rownames = FALSE))
    }

    if(lang == 'html') rendHTML(R) else rendHTML(R, html = FALSE)
  }

## html.validate: thin delegate, matching html.anova.rms/html.summary.rms's
## identical pattern -- lang detection inside latex.validate (via
## prType()) already handles the branching.
html.validate <- function(object, digits=4, B=Inf, caption=NULL, ...)
  latex.validate(object, digits=digits, B=B, caption=caption, ...)

## -----------------------------------------------------------------------
## Detailed notes on extending validate.ols's table output to typst,
## referenced by the short [N] tags inline above. Most of this is
## simply building out Typst coverage for the first time -- part of
## this project's goal from the start, not corrections to anything
## broken. Two items are worth reading carefully regardless (notes [1]
## and [6]): one is a bug class that has bitten this exact spot before
## in a sibling file, the other is a real, deliberate simplification of
## existing parameters.
##
## [1] print.validate previously only auto-dispatched for html
##     (if(prType()=='html') return(html.validate(...))); latex and
##     typst fell through to plain-text printing entirely, even under
##     prType='latex'. Extended to match html's behavior, per explicit
##     confirmation this was wanted (unlike anova.rms/summary.rms, where
##     the equivalent auto-dispatch already existed for all non-plain
##     languages). Every branch uses return() -- the exact same bug
##     class already found and fixed once in summary.rms.s: rendHTML()'s
##     asis_output mechanism requires being the visible, top-level
##     return value to trigger at all, and a branch without return()
##     would have its output silently discarded when the function
##     continued on past the switch/if-chain, producing no visible
##     output despite no error. Worth remembering for any future
##     function converted to this mechanism.
##
## [2] latex.validate/html.validate unification: previously two entirely
##     separate, fully duplicated functions -- one manually writing
##     \Needspace{}/\begin{center}/Hmisc::latex()/\end{center} blocks
##     directly to a file connection via cat(), the other building an R
##     character vector via htmlTable::htmlTable() calls and emitting
##     once via rendHTML(). Unified into one lang-aware implementation
##     using rms_tiny_table() (shared with anova.rms and summary.rms,
##     extended for this file's needs -- see anova.rms.s's own note [6]
##     for the show_rownames/linebreak additions), with html.validate
##     reduced to a thin delegate, matching the pattern already
##     established for html.anova.rms/html.summary.rms.
##
## [3] Row-name translation (Dxy, R2, Emax, g, gp, gr, rho, pdm, etc.):
##     extended with a third, typst-specific translation vector,
##     following the same conventions already established and
##     compile-tested in prStats' trans_tbl -- quoted multi-letter
##     subscripts ($D_("xy")$), bare greek names ($rho$), bare
##     single-letter symbols left unquoted. The pdm entry (the mean
##     |P(Y>=median)-0.5| statistic) is the exact same underlying
##     statistic as prStats' '|Pr(Y>=median)-0.5|' row, so the identical
##     typst expression already built there is reused here rather than
##     rebuilt independently -- still the single highest-risk expression
##     in this whole effort (a long expression inside overline(), never
##     itself directly compile-tested), worth testing deliberately here
##     too, not just assumed correct because it matches prStats'.
##
##     Emax's subscript ("max") is quoted ($E_("max")$) rather than left
##     bare, even though "max" is a recognized Typst function/operator
##     name -- its behavior as bare subscript content (not a function
##     call) is genuinely uncertain, and quoting guarantees it displays
##     as literal text regardless.
##
## [4] Bullet marker for retained variables (previously '$\\bullet$'
##     for latex, htmlSpecial('mediumsmallwhitecircle') for html): typst
##     uses the literal Unicode bullet character (\u2022) directly,
##     following the same "prefer literal Unicode over a guessed Typst
##     symbol name" convention already used elsewhere in this project
##     (e.g. the >=/<= comparison operators in prStats).
##
## [5] Vestigial parameters, kept for signature compatibility with
##     existing direct calls but no longer used internally: file/append
##     (the file-writing mechanism is replaced entirely by rendHTML()'s
##     asis-output mechanism, matching every other function converted
##     this way in this project); title (only ever used by the old
##     Hmisc::latex() calls, which rms_tiny_table() has no equivalent
##     parameter for -- also worth noting the original default,
##     first.word(deparse(substitute(x))), referenced a parameter named
##     x that doesn't exist in this function's own signature (it's
##     named object), a pre-existing bug unrelated to typst; harmless
##     here specifically because R's lazy evaluation means this default
##     expression is never evaluated now that title is unused); size/
##     extracolsize (LaTeX-specific font-sizing directives with no
##     equivalent needed in the tinytable-based approach).
##
## [6] table.env's former role (choosing between Hmisc::latex()'s own
##     numbered-caption mechanism and a manually cat()'d plain-text
##     caption) is no longer used -- captions are always shown via
##     rms_tiny_table()'s plain bold-heading mechanism when given,
##     matching the simplification already applied in anova.rms.s and
##     summary.rms.s. This is a direct match to the ORIGINAL code's net
##     behavior, not a new simplification specific to this file: the
##     original latex.validate showed the caption either way (via
##     Hmisc::latex()'s caption= when table.env=TRUE, or via a manual
##     cat() when table.env=FALSE) -- table.env only ever chose the
##     mechanism, never whether the caption appeared at all.
## -----------------------------------------------------------------------
