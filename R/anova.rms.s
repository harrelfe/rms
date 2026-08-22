#main.effect=F to suppress printing main effects when the factor in
#question is involved in any interaction.

anova.rms <- function(object, ..., main.effect=FALSE, tol=.Machine$double.eps,
                      test=c('F', 'Chisq', 'LR'),
                      india=TRUE, indnl=TRUE, ss=TRUE,
                      vnames=c('names', 'labels'),
                      posterior.summary=c('mean', 'median', 'mode'),
                      ns=500, cint=0.95, fitargs=NULL) {

  misstest   <- missing(test)
  test       <- match.arg(test)
  fitfunName <- class(object)


  if(test == 'LR') {
    if('x' %nin% names(object))
      stop('must specify x=TRUE, y=TRUE when fitting model to use LR test')
    X <- object[['x']]
    fitter <- do.call('quickRefit', c(list(object, what='fitter'), fitargs))
    # Compute deviance for full model
    devf <- getDeviance(object)
    devf <- devf[length(devf)]
  }

 ava <- if(test == 'LR') function(idx) {
   if(length(idx) == ncol(X)) c(object$stats['Model L.R.'], length(idx))
   else {
    devs <- getDeviance(fitter(X[, - idx, drop=FALSE]), fitfunName)
    devs <- devs[length(devs)]
    c(devs - devf, length(idx))
   }
  } else function(idx) {
    chisq <- coef[idx] %*% solve(cov[idx, idx], coef[idx], tol=tol)
    c(chisq, length(idx))
  }
             

  eEV <- function(test=integer()) {
    coef  <- if(length(test)) draws[, test, drop=FALSE] else draws
    co    <- if(length(test)) cov[test, test, drop=FALSE] else cov
    m     <- nrow(coef)
    chisq <- numeric(m)
    for(i in 1 : m)
      chisq[i] <- coef[i,, drop=FALSE] %*%
        solvet(co, t(coef[i,, drop=FALSE]), tol=tol)

    if(! length(test)) return(chisq)   # Assumes stored in chisqT

    ## variance explained by a variable/set of variables is
    ## approximated by the Wald chi-square
    ## pev = partial explained variation = chisq/(chisq for full model)
    pev <- chisq / chisqT
    ## Overall pev is the pev at the posterior mean/median beta (last element)
    ## Also compute HPD interval.
    ci <- rmsb::HPDint(pev[-m], cint)
    c(REV=pev[m], Lower=ci[1], Upper=ci[2], d.f.=length(test))
  }

  obj.name <- as.character(sys.call())[2]
  vnames <- match.arg(vnames)
  posterior.summary <- match.arg(posterior.summary)
  is.ols <- inherits(object,'ols')

  nrp <- num.intercepts(object)
  cov <- vcov(object, regcoef.only=TRUE, intercepts='none')

  draws <- object$draws
  bayes <- length(draws) > 0
  chisqBayes <- NULL

  if(bayes) {
    if(nrp > 0) draws <- draws[, -(1 : nrp), drop=FALSE]

    betaSummary <- rmsb::getParamCoef(object, posterior.summary)
    if(nrp > 0) betaSummary <- betaSummary[-(1 : nrp)]

    X <- object[['x']]
    if(! length(X)) stop('x=TRUE must have been specified to fit')
    nc <- ncol(X)
    ndraws <- nrow(draws)
    ns     <- min(ndraws, ns)
    if(ns < ndraws) {
      j <- sample(1 : ndraws, ns, replace=FALSE)
      draws <- draws[j,, drop=FALSE]
    }

    ## Augment draws with a last row with posterior central tendency
    draws <- rbind(draws, posteriorSummary=betaSummary)
    ## Compute variances of linear predictor without omitting variables
    chisqT <- eEV()
    m  <- length(chisqT)
    ci <- rmsb::HPDint(chisqT[-m], cint)
    chisqBayes <- c(chisqT[m], ci)
    names(chisqBayes) <- c('Central', 'Lower', 'Upper')
  }

  if(misstest) test <- if(is.ols) 'F' else 'Chisq'
  if(! is.ols && test=='F') stop('F-test not allowed for this type of model')
  if(bayes) test <- 'Chisq'

  if(! is.ols) ss <- FALSE

  at     <- object$Design
  assign <- object$assign
  name   <- at$name
  labels <- at$label
  nama   <- names(assign)[1]
  asso   <- 1*(nama=="(Intercept)" | nama=="Intercept")
  names(assign)[-asso] <- name
  namelab <- if(vnames == 'names') name else labels

  ia <- at$interactions
  nia <- if(!length(ia)) 0 else ncol(ia)

  assume <- at$assume.code
  parms  <- at$parms
  f      <- length(assume)

  ## If using labels instead of names, substitute labels in interaction labels,
  ## e.g. change x1 * x2 to label(x1) * label(x2)
  if(vnames == 'labels' && any(assume == 9)) {
    for(i in which(assume == 9)) {
      parmi <- parms[[name[i]]]
      parmi <- parmi[, 1][parmi[, 1] > 0]
      namelab[i] <- paste(labels[parmi], collapse=' * ')
    }
  }

  ncall <- names(sys.call())[-(1 : 2)]
  alist <- as.character(sys.call())[-(1 : 2)]
  if(length(alist) && length(ncall)) alist <- alist[ncall == '']
  which <- if(length(alist)) {
    jw <- charmatch(alist, name, 0)
    if(any(jw == 0))
      stop(paste("factor names not in design: ",
                 paste(alist[jw == 0], collapse=" ")))
    jw
  }
  else 1 : f

  if(! bayes) {
    if(length(object$est) && !length(object$u))
      stop("est in fit indicates score statistics but no u in fit")

    if(test != 'LR' && ! length(object$coefficients))
        stop("estimates not available for Wald statistics")

    coef <- object$coefficients
    cik  <- attr(coef, 'intercepts')
#    }
#    else {
#      if(!length(object$u))
#        stop("score statistics not available")
#      coef <- object$u
#    }
  }

cov <- vcov(object, regcoef.only=TRUE, intercepts='none')

  if(bayes) for(j in 1:length(assign))
        assign[[j]] <- assign[[j]] - nrp

  else {
    ## Omit row/col for scale parameters
    ## Compute # intercepts nrpcoef to skip in testing
    nrpcoef <- num.intercepts(object, 'coef')
    if(nrpcoef > 0) {
      coef <- coef[-(1 : nrpcoef)]
      for(j in 1:length(assign))
        assign[[j]] <- assign[[j]] - nrpcoef
    }

    nc <- length(coef)
  }

  dos <- if(bayes) eEV else ava

  stats <- NULL
  lab   <- NULL
  W     <- vinfo <- list()
  s     <- 0
  all.slopes <- rep(FALSE, nc)
  all.ia     <- rep(FALSE, nc)
  all.nonlin <- rep(FALSE, nc)
  num.ia     <- 0
  num.nonlin <- 0
  issue.warn <- FALSE

  for(i in which) {
    j <- assume[i]
    parmi <- parms[[name[i]]]
    low.fact <- if(j != 9) i else (parmi[,1])[parmi[,1] > 0]

    nl <- if(!length(names(at$nonlinear))) at$nonlinear[[i]]
    else at$nonlinear[[name[i]]]

    if(!length(nl))
      nl <- rep(FALSE, length(assign[[name[i]]]))

    ## Factor no. according to model matrix is 1 + number of non-strata factors
    ## before this factor
    if(j != 8) {
      ##ignore strata
      jfact <- if(i==1) 1 else 1 + sum(assume[1 : (i - 1)] != 8)

      main.index <- assign[[jfact + asso]]
      nonlin.ia.index <- NULL	#Should not have to be here. Bug in S?
      all.slopes[main.index] <- TRUE
      ni <- if(nia == 0) 0 else sum(ia == i)

      if(nia==0) ni <- 0
      else
        for(k in 1:ncol(ia))
          ni <- ni + !any(is.na(match(low.fact, ia[,k])))

      if(ni==0 | main.effect) {
        w <- dos(main.index)
        s <- s+1; W[[s]] <- main.index
        stats <- rbind(stats, w)
        lab <- c(lab, namelab[i])
        vinfo[[s]] <- list(name=name[i], type='main effect')
      }

      ## If term is involved in any higher order effect, get pooled test
      ## by adding in all high-order effects containing this term
      ## For 2nd order interaction, look for 3rd order interactions
      ## containing both factors
      ## nonlin.ia.index <- NULL	#Used to be here.  Bug in S?
      if(ni > 0) {
        ia.index <- NULL
        mm <- (1:f)[assume == 9]
        mm <- mm[mm != i]
        for(k in mm) {
          parmk <- parms[[name[k]]]
          hi.fact <- parmk[,1]
          m <- match(low.fact, hi.fact)
          if(!any(is.na(m))) {
            kfact <- if(k==1) 1 else
            1 + sum(assume[1:(k-1)] != 8)

            idx <- assign[[kfact + asso]]
            ia.index <- c(ia.index, idx)

            if(ncol(parmk)>1)
              for(jj in 1:length(m))
                nonlin.ia.index <- c(nonlin.ia.index,
                                     idx[parmk[m[jj],-1] == 1])

            nonlin.ia.index <- if(length(nonlin.ia.index))
              unique(nonlin.ia.index)
            else NULL
            ##Highest order can be counted twice
          }
        }

        idx <- c(main.index, ia.index)
        all.slopes[idx] <- TRUE
        w <- dos(idx)
        s <- s + 1; W[[s]] <- idx
        stats <- rbind(stats, w)
        lab <- c(lab, paste(namelab[i],
                            " (Factor+Higher Order Factors)"))
        vinfo[[s]] <- list(name=name[low.fact],
                           type=if(j==9) 'interaction' else 'combined effect')

        ## If factor i in >1 interaction, print summary
        ## Otherwise, will be printed later
        if(india && (j != 9 & ni > 1)) {
          w <- dos(ia.index)
          s <- s + 1; W[[s]] <- ia.index
          stats <- rbind(stats, w)
          lab <- c(lab, " All Interactions")
          vinfo[[s]] <- list(name=name[low.fact], type='combined interactions')
        }
      }
      if(any(nl) && any(!nl)) {
        ## Tests of adequacy of linear relationship
        idx <- c(main.index[nl], nonlin.ia.index)

        num.nonlin <- num.nonlin+1
        all.nonlin[idx] <- TRUE
        if(indnl) {
          w <- dos(idx)
          s <- s + 1; W[[s]] <- idx
          stats <- rbind(stats, w)
          lab <- c(lab, if(!length(nonlin.ia.index))" Nonlinear"
          else " Nonlinear (Factor+Higher Order Factors)")
          vinfo[[s]] <- list(name=name[low.fact],
                             type=if(j==9)
                             'nonlinear interaction' else 'nonlinear')
        }
      }
      ## If interaction factor involves a non-linear term from an
      ## expanded polynomial, lspline, rcspline, or scored factor,
      ## do tests to see if a simplification (linear interaction) is
      ## adequate.  Do for second order only.
      if(j == 9) {
        num.ia <- num.ia+1
        all.ia[main.index] <- TRUE
        if(parmi[3,1] > 0)
          issue.warn <- TRUE

        if(parmi[3,1] == 0 && ncol(parmi) > 1) {
          nonlin.x <- as.logical(parmi[1,2:ncol(parmi)])
          nonlin.y <- as.logical(parmi[2,2:ncol(parmi)])
          nonlin.xy <- nonlin.x | nonlin.y
          nonlin.xandy <- nonlin.x & nonlin.y
          idx <- main.index[nonlin.xy]
          li <- length(idx)

          if(li > 0) {
            num.nonlin <- num.nonlin+1
            all.nonlin[idx] <- TRUE
            if(indnl) {
              w <- dos(idx)
              s <- s + 1
              W[[s]] <- idx
              stats  <- rbind(stats, w)
              lab    <- c(lab," Nonlinear Interaction : f(A,B) vs. AB")
              vinfo[[s]] <- list(name=name[low.fact],
                                 type='nonlinear interaction')
            }
            idx <- main.index[nonlin.xandy]
            li <- length(idx)

            if(indnl && li > 0) {
              w <- dos(idx)
              s <- s + 1
              W[[s]] <- idx
              stats <- rbind(stats,w)
              lab <- c(lab, " f(A,B) vs. Af(B) + Bg(A)")
              vinfo[[s]] <- list(name=name[low.fact],
                                 type='doubly nonlinear interaction')
            }

            idx <- main.index[nonlin.x]
            li <- length(idx)
            if(indnl && (li > 0 & any(nonlin.x != nonlin.xy))) {
              w <- dos(idx)
              s <- s+1
              W[[s]] <- idx
              stats <- rbind(stats, w)
              lab   <- c(lab, paste(" Nonlinear Interaction in",
                                    namelab[parmi[1,1]],"vs. Af(B)"))
              vinfo[[s]] <-
                list(name=name[low.fact],
                     type='nonlinear interaction in first variable')
            }

            idx <- main.index[nonlin.y]
            li <- length(idx)

            if(indnl && (li > 0 & any(nonlin.y != nonlin.xy))) {
              w <- dos(idx)
              s <- s + 1
              W[[s]] <- idx
              stats <- rbind(stats,w)
              lab <- c(lab, paste(" Nonlinear Interaction in",
                                namelab[parmi[2,1]],"vs. Bg(A)"))
              vinfo[[s]] <-
                list(name=name[low.fact],
                     type='nonlinear interaction in second variable')
            }
          }
        }
      }
    }
  }
  ## If all lines so far had (Factor +Higher Order Factors) in them,
  ## remove this redundancy
  if(length(grep('\\(Factor\\+Higher Order Factors\\)', lab)) == length(lab))
    lab <- gsub('\\(Factor\\+Higher Order Factors\\)', '', lab)

  ## If >1 test of adequacy, print pooled test of all nonlinear effects
  if(num.nonlin > 1 || (num.nonlin==1 & !indnl)) {
    idx <- (1:nc)[all.nonlin]
    li <- length(idx)
    w <- dos(idx)
    s <- s + 1; W[[s]] <- idx
    stats <- rbind(stats, w)
    lab <- c(lab, "TOTAL NONLINEAR")
    vinfo[[s]] <- list(type='total nonlinear')
  }

  ## If >1 test of interaction, print pooled test of all interactions in list
  if(num.ia > 1 || (num.ia==1 & !india)) {
    idx <- (1:nc)[all.ia]
    li <- length(idx)
    w <- dos(idx)
    s <- s+1
    W[[s]] <- idx
    stats <- rbind(stats, w)
    lab <- c(lab, "TOTAL INTERACTION")
    vinfo[[s]] <- list(type='total interaction')
  }

  ## If >0 test of adequacy and >0 test of interaction, print pooled test of
  ## all nonlinear and interaction terms
  if(num.nonlin > 0 & num.ia > 0) {
    idx <- (1:nc)[all.nonlin | all.ia]
    li <- length(idx)
    w <- dos(idx)
    s <- s + 1
    W[[s]] <- idx
    stats <- rbind(stats,w)
    lab <- c(lab, "TOTAL NONLINEAR + INTERACTION")
    vinfo[[s]] <- list(type='complexity')
  }

  ## Get total test for all factors listed
  idx <- (1:nc)[all.slopes | all.ia]
  w   <- dos(idx)
  s   <- s + 1;  W[[s]] <- idx
  stats <- rbind(stats, w)
  lab <- c(lab, "TOTAL")
  vinfo[[s]] <- list(type='global')

statnam <- if(bayes) c('REV', 'Lower', 'Upper', 'd.f.')
           else
             c('Chi-Square','d.f.')

  if(! bayes) {
    if(is.ols) {
      sigma2 <- object$stats['Sigma']^2
      dfe    <- object$df.residual
    }

    if(ss) {
      stats <- cbind(stats[,2], stats[,1]*sigma2, stats[,1]*sigma2/stats[,2],
                     stats[,1])
      statnam <- c('d.f.', 'Partial SS', 'MS', 'Chi-Square')
      stats <- rbind(stats, Error=c(dfe, sigma2*dfe, sigma2, NA))
      s <- s + 1; W[[s]] <- NA
      lab <- c(lab, 'ERROR')
      vinfo[[s]] <- list(type='error')
    }

    j <- statnam == 'Chi-Square'
    dfreg <- stats[, statnam=='d.f.']

    if(test == 'F') {
      stats[,j] <- stats[,j] / dfreg
      statnam[j] <- 'F'
      stats <- cbind(stats, P=1 - pf(stats[,j], dfreg, dfe))
      attr(stats,'df.residual') <- dfe
    }
    else
      stats <- cbind(stats,1 - pchisq(stats[,j], dfreg))

    statnam <- c(statnam, 'P')
  }

  dimnames(stats) <- list(lab, statnam)
  attr(stats,'formula') <- formula(object)
  attr(stats,"obj.name") <- obj.name
  attr(stats,"class") <- c("anova.rms","matrix")

  names(W) <- lab
  attr(stats, 'which') <- W
  attr(stats, 'test')  <- test
  if(! bayes) attr(stats,"coef.names") <- names(coef)
  attr(stats,"non.slopes") <- nrp
  attr(stats,"vinfo")      <- vinfo
  attr(stats,"chisqBayes") <- chisqBayes

  if(issue.warn)
    warning("tests of nonlinear interaction with respect to single component \nvariables ignore 3-way interactions")

  stats
}

print.anova.rms <- function(x, which=c('none','subscripts',
                                       'names','dots'),
                            table.env=FALSE,
                            ...) {
  which <- match.arg(which)
  lang  <- prType()

  stats <- x
  digits <- c('Chi-Square'=2, F=2, 'd.f.'=0, 'Partial SS'=15, MS=15, P=4,
              REV=3, Lower=3, Upper=3)
  cstats <- matrix('', nrow=nrow(stats), ncol=ncol(stats),
                   dimnames=dimnames(stats))

  bchi <- attr(stats, 'chisqBayes')

  test <- attr(stats, 'test')
  if(! length(test)) test <- 'Chisq'   # for legacy fit objects
  if(test == 'LR')    test <- 'Likelihood Ratio'
  if(test == 'Chisq') test <- 'Wald'

  do.which <- which!='none' && length(W <- attr(stats,'which'))
  params <- NULL

  if(do.which) {
    if(which=='subscripts')
      simplifyr <- function(x) {
        x <- sort(unique(x))
        n <- length(x)
        ranges <- character(n)
        m <- 0
        s <- x

        while(length(s) > 0) {
          j <- s == s[1] + (1:length(s))-1
          m <- m+1
          ranges[m] <- if(sum(j)>1) paste(range(s[j]),collapse='-')
          else s[1]
          s <- s[!j]
        }

        ranges[1:m]
      }

    k <- length(W)
    w <- character(k)
    coef.names <- attr(stats,'coef.names')
    for(i in 1:k) {
      z <- W[[i]]

      if(all(is.na(z))) w[i] <- ''
      else {
        z <- sort(z)
        w[i] <- switch(which,
                       subscripts=paste(simplifyr(z), collapse=','),
                       names=paste(coef.names[z],collapse=','),
                       dots={
                         dots <- rep(' ',length(coef.names))
                         dots[z] <- '.'
                         paste(dots, collapse='')
                       })
      }
    }
    params <- w
    if(lang == 'html' && which == 'dots') {
      params <- gsub(' ',   '&nbsp;', params)     # non-breaking space
      params <- gsub('\\.', '\u2022', params)     # bullet
      params <- paste0('<tt>', params, '</tt>')   # monospace
      }
  }   # end do.which

  if(lang != 'plain')
    return(latex.anova.rms(x, file='', table.env=table.env,
                           params=params, ...))

  sn <- colnames(cstats)

  for(j in 1:ncol(cstats))
    cstats[,j] <- format(round(stats[,j], digits[sn[j]]))

  cstats[is.na(stats)] <- ''
  j <- sn=='P'
  cstats[stats[,j] < 0.00005,j] <- '<.0001'
  cstats <- cbind(rownames(stats), cstats)

  dimnames(cstats) <- list(rep("",nrow(stats)),
                           c("Factor    ",colnames(stats)))

  heading <- if(length(bchi))
        paste('                  Relative Explained Variation        Response:',
              as.character(attr(stats, "formula")[2]), sep = "")
             else
             paste("                ",
                   if(any(colnames(stats) == 'F')) "Analysis of Variance"
                   else paste(test, "Statistics"),
                   "          Response: ",
                   as.character(attr(stats, "formula")[2]), sep = "")

  cat(heading,"\n\n")
  if(any(sn=='MS'))
    cstats[cstats[,1]=='TOTAL',1] <- 'REGRESSION'

  if(do.which) cstats <- cbind(cstats, Tested=w)

  print(cstats, quote=FALSE)

  if(do.which && which!='names') {
    cat('\nSubscripts correspond to:\n')
    print(coef.names, quote=FALSE)
  }

  if(!any(sn=='MS') && length(dfe <- attr(stats,'df.residual')))
    cat('\nError d.f.:', dfe, '\n')

  if(length(bchi)) {
    bchi <- round(bchi, 1)
    cat('\nApproximate total model Wald total chi-square used in denominators of REV:\n',
        bchi['Central'], ' [', bchi['Lower'], ', ',
        bchi['Upper'], ']\n', sep='')
  }

  invisible()
}

## rms_tiny_table: renders d (a data.frame whose rownames are the
## desired row labels) as a table for the current output language, with
## an optional bold centered caption above it. Used by anova.rms,
## summary.rms, and validate.ols. Row labels become an explicit
## left-justified first column; data columns right-justified; header
## centered; rows striped.
##
## Structured as one format-agnostic core builder plus a separate,
## named finalize/caption function per language, so each format's
## handling is self-contained and can be extended independently. See
## end of file for the detailed reasoning behind each format's specific
## choices.
## show_rownames=FALSE: omits the row-label column entirely, for tables
## with no meaningful row names (validate.ols's marker/frequency
## tables). linebreak="\n" (unconditional): multi-line column headers,
## harmless for headers without embedded newlines. align_data/
## bold_header/header_rule/vertical_borders/width: new, all defaulting
## to prior behavior exactly (right-aligned data, plain non-bold header,
## no rule, no borders, no explicit width) -- added to absorb
## prStats/stats_tiny_table's styling needs, previously a separate,
## partially-duplicated implementation. See note [6] at end of file.
rms_tiny_table_core <- function(d, lang, fontsize = NULL,
                                show_rownames = TRUE, align_data = 'r',
                                bold_header = FALSE, header_rule = FALSE,
                                vertical_borders = FALSE, width = NULL) {
  dd <- if(show_rownames) {
    rowlab <- rownames(d)
    tmp <- data.frame(Var = rowlab, unclass(d), check.names = FALSE)
    colnames(tmp) <- c('', colnames(d))
    tmp
  } else if(is.data.frame(d)) d
    else as.data.frame(unclass(d), check.names = FALSE)
  ## is.data.frame(d) check added -- see end-of-file note on the
  ## stats_tiny_table regression this fixed: the unclass()/as.data.frame()
  ## round-trip was new (the original stats_tiny_table called tt()
  ## directly on its data.frame, no conversion step), and produced
  ## corrupted output specifically for stats_tiny_table's multi-line,
  ## marked-up column headers -- confirmed from a real compile error.
  ## Every current caller already passes a proper data.frame here, so
  ## skipping the round-trip when one is already given sidesteps the
  ## issue without depending on fully diagnosing the exact interaction.
  nc <- ncol(dd)
  nr <- nrow(dd)
  jdata <- if(show_rownames) 2:nc else 1:nc

  x <- if(length(width)) tinytable::tt(dd, output = lang, width = width)
       else tinytable::tt(dd, output = lang)
  x <- tinytable::format_tt(x, j = 1:nc, escape = FALSE)
  x <- tinytable::format_tt(x, i = 0,    escape = FALSE, linebreak = "\n")
  if(show_rownames) x <- tinytable::style_tt(x, j = 1, align = 'l')
  ## align_data as a vector: applied one column at a time rather than
  ## passing the vector directly to a single style_tt(j=jdata, align=)
  ## call, since element-wise recycling of align across j hasn't been
  ## confirmed for tinytable's style_tt() -- looping is the safe,
  ## already-proven pattern used everywhere else in this file (e.g. the
  ## j=1/jdata split just above). A scalar align_data (the prior,
  ## default behavior) is unaffected -- single style_tt() call, same as
  ## before.
  if(length(align_data) > 1)
    for(jj in seq_along(jdata))
      x <- tinytable::style_tt(x, j = jdata[jj], align = align_data[jj])
  else
    x <- tinytable::style_tt(x, j = jdata, align = align_data)
  x <- tinytable::style_tt(x, i = 0,     align = 'c', bold = bold_header)
  if(header_rule) x <- tinytable::style_tt(x, i = 0, j = 1:nc, line = "b")
  if(vertical_borders) {
    x <- tinytable::style_tt(x, i = 0:nr, j = 1:nc, line = "l")
    x <- tinytable::style_tt(x, i = 0:nr, j = nc,   line = "r")
  }
  x <- tinytable::theme_striped(x)
  if(length(fontsize) && fontsize != 1)
    x <- tinytable::style_tt(x, fontsize = fontsize)
  x
}

## inner=: optional raw tabularray inner options (e.g. colsep/rowsep
## tuning), new, defaulting to NULL (prior behavior: no inner= passed).
rms_tiny_table_finalize_latex <- function(x, inner = NULL) {
  x <- if(length(inner))
         tinytable::theme_latex(x, environment_table = FALSE, inner = inner)
       else
         tinytable::theme_latex(x, environment_table = FALSE)
  result <- tinytable::save_tt(x, 'latex')
  paste0('\\begin{center}\n', result, '\n\\end{center}')
}

rms_tiny_table_finalize_html <- function(x) tinytable::save_tt(x, 'html')

rms_tiny_table_finalize_typst <- function(x) {
  result <- tinytable::save_tt(x, 'typst')
  typstFence <- function(s) paste0("````{=typst}\n", s, "\n````\n")
  typstFence(paste0('#block(width: 100%)[#align(center)[\n', result, '\n]]'))
}

rms_tiny_table_caption_latex <- function(caption)
  paste0('\\begin{center}\\textbf{', caption, '}\\end{center}\n')

rms_tiny_table_caption_html <- function(caption)
  paste0('<p style="text-align:center"><strong>', caption, '</strong></p>\n')

rms_tiny_table_caption_typst <- function(caption) {
  typstFence <- function(s) paste0("````{=typst}\n", s, "\n````\n")
  typstFence(paste0('#block(width: 100%)[#align(center)[#strong[',
                    caption, ']]]'))
}

rms_tiny_table <- function(d, lang, caption = NULL, fontsize = NULL,
                           show_rownames = TRUE, align_data = 'r',
                           bold_header = FALSE, header_rule = FALSE,
                           vertical_borders = FALSE, width = NULL,
                           latex_inner = NULL) {
  x <- rms_tiny_table_core(d, lang, fontsize = fontsize,
                           show_rownames = show_rownames,
                           align_data = align_data,
                           bold_header = bold_header,
                           header_rule = header_rule,
                           vertical_borders = vertical_borders,
                           width = width)
  result <- switch(lang,
                   latex = rms_tiny_table_finalize_latex(x, inner = latex_inner),
                   html  = rms_tiny_table_finalize_html(x),
                   typst = rms_tiny_table_finalize_typst(x))
  if(length(caption)) {
    capMarkup <- switch(lang,
                        latex = rms_tiny_table_caption_latex(caption),
                        html  = rms_tiny_table_caption_html(caption),
                        typst = rms_tiny_table_caption_typst(caption))
    result <- c(capMarkup, result)
  }
  result
}

latex.anova.rms <-
  function(object,
           title=paste('anova', attr(object, 'obj.name'), sep='.'),
           dec.chisq=2, dec.F=2, dec.ss=NA,
           dec.ms=NA, dec.P=4, dec.REV=3,
           table.env=TRUE, caption=NULL,
           fontsize=1, params=NULL, ...) {

    ## params is only used if called from print.anova.rms
    ## It is not meant to be provided by the user in a latex. call

    lang <- prType()
    html <- lang == 'html'

    sn   <- colnames(object)
    rowl <- rownames(object)
    if(any(sn=='MS')) rowl[rowl=='TOTAL'] <- 'REGRESSION'

    ## extended to typst -- see note [1] at end of file
    rowl <- switch(lang, latex = latexTranslate(rowl),
                   typst = typstTranslate(rowl), rowl)

    specs <- markupSpecs[[lang]]
    bold  <- specs$bold
    math  <- specs$math

    ## Translate interaction symbol (*) to times symbol
    ## rowl <- gsub('\\*', specs$times, rowl)   # changed * to $times$
    ## extended to typst -- see note [2] at end of file
    pat <- if(lang == 'typst') '\\*' else '*'
    rowl <- gsub(pat, specs$times, rowl, fixed=TRUE)

    ## Put TOTAL rows in boldface
    rowl <- ifelse(substring(rowl, 1, 5) %in% c("REGRE", "ERROR", "TOTAL"),
                   bold(rowl), rowl)

    rowl <- ifelse(substring(rowl, 1, 1) == " ",
                 paste0(specs$lspace, specs$italics(substring(rowl,2)), sep=""),
                 rowl) # preserve leading blank

    P <- object[,3]

    dstats <- as.data.frame(object)
    attr(dstats, 'row.names') <- rowl

    digits <- c('Chi-Square'=dec.chisq, F=dec.F, 'd.f.'=0,
                'Partial SS'=dec.ss, MS=dec.ms, P=dec.P,
                REV=dec.REV, Lower=dec.REV, Upper=dec.REV)

    dig <- digits[sn]
    ## column names double-processed for the chi-square entry -- see
    ## note [4] at end of file. specs$chisq()'s own $-wrapping is
    ## format-inconsistent (see note [5]), checked directly below rather
    ## than assumed.
    nms <- ifelse(sn %nin% c('d.f.','MS','Partial SS','Chi-Square'),
                  math(sn), sn)
    chi2lab <- specs$chisq(add='')
    nms[sn == 'Chi-Square'] <-
      if(substring(chi2lab, 1, 1) == '$') chi2lab else math(chi2lab)
    names(dstats) <- nms

    resp <- as.character(attr(object, 'formula')[2])
    ## extended to typst -- see note [1] at end of file
    resp <- switch(lang, latex = latexTranslate(resp),
                    typst = typstTranslate(resp), resp)

    test <- attr(object, 'test')
    if(! length(test))  test <- 'Chisq'   # for legacy fit objects
    if(test == 'LR')    test <- 'Likelihood Ratio'
    if(test == 'Chisq') test <- 'Wald'

    bchi <- attr(object, 'chisqBayes')
    wl <- if(length(bchi)) 'Relative Explained Variation' else
      if(any(sn == 'F')) 'Analysis of Variance' else paste(test, 'Statistics')
    if(! length(caption))
      caption <- paste0(wl, " for ", specs$code(resp))

    i <- 0
    for(nn in names(dstats)) {
      i <- i + 1
      dstats[[nn]] <- formatNP(dstats[[nn]], digits=dig[i],
                               lang   = lang,
                               pvalue = nn == math('P'))
    }

  if(length(bchi)) {
    bchi <- round(bchi, 1)
    w <- paste0('Approximate total model Wald ',
                specs$math(specs$chisq(add='')),
               ' used in denominators of REV:',
        bchi['Central'], ' [', bchi['Lower'], ', ',
        bchi['Upper'], '].')
    caption <- paste0(caption, '. ', w)
  }

    if(length(params)) {
      dstats$Tested <- params
      sn <- c(sn, 'Tested')
    }

    ## table-making via tinytable, modular across latex/html/typst --
    ## see rms_tiny_table() above; see note [3] at end of file
    w <- rms_tiny_table(dstats, lang, caption = caption, fontsize = fontsize)
    if(lang == 'html') rendHTML(w) else rendHTML(w, html = FALSE)
  }


html.anova.rms <-
  function(object, ...) latex.anova.rms(object, ...)


plot.anova.rms <-
  function(x, what=c("chisqminusdf","chisq","aic",
                "P","partial R2","remaining R2",
                "proportion R2", "proportion chisq"),
           xlab=NULL,
           pch=16, rm.totals=TRUE, rm.ia=FALSE,
           rm.other=NULL, newnames,
           sort=c("descending","ascending","none"),
           margin=c('chisq', 'P'),
           pl=TRUE, trans=NULL, ntrans=40, height=NULL, width=NULL, ...) {

    what <- match.arg(what)
    sort <- match.arg(sort)
    isbase <- Hmisc::grType() == 'base'


    htmlSpecs <- markupSpecs$html
    schisq    <- htmlSpecs$chisq()
    nbsp      <- htmlSpecial('nbsp')

    if(! length(xlab)) {

      xlab <-
        if(isbase)
          switch(what,
                 chisq=expression(chi^2),
                 "proportion chisq"=expression(paste("Proportion of Overall", ~chi^2)),
                 chisqminusdf=expression(chi^2~-~df),
                 aic="Akaike Information Criterion",
                 P="P-value",
                 "partial R2"=expression(paste("Partial",~R^2)),
                 "remaining R2"=expression(paste("Remaining~",R^2,
                                                 "~After Removing Variable")),
                 "proportion R2"=expression(paste("Proportion of Overall",
                                                  ~R^2)))
        else
          switch(what,
                 chisq        = schisq,
                 "proportion chisq" = paste('Proportion of Overall', schisq),
                 chisqminusdf = paste0(schisq, nbsp, '-', nbsp, 'df'),
                 aic          = "Akaike Information Criterion",
                 P            = "P-value",
                 "partial R2" = 'Partial R<sup>2</sup>',
                 "remaining R2" = 'Remaining R<sup>2</sup> After Removing Variable',
                 "proportion R2"='Proportion of Overall R<sup>2</sup>')
      }
    rm <- c(if(rm.totals) c("TOTAL NONLINEAR","TOTAL NONLINEAR + INTERACTION",
                            "TOTAL INTERACTION","TOTAL"),
            " Nonlinear"," All Interactions", "ERROR",
            " f(A,B) vs. Af(B) + Bg(A)", rm.other)

    rn <- rownames(x)
    rm <- c(rm, rn[substring(rn, 2, 10) == "Nonlinear"])
    k <- !(rn %in% rm)
    if(rm.ia) k[grep("\\*", rn)] <- FALSE

    an <- x[k,, drop=FALSE]

    if(! isbase && ! length(height))
      height <- plotlyParm$heightDotchart(nrow(an))

    if('REV' %in% colnames(x)) {    # Bayesian
      xlab <- 'Relative Explained Variation'
      i <- switch(sort,
                  none       = 1 : nrow(an),
                  descending = order(an[, 'REV'], decreasing=TRUE),
                  ascending  = order(an[, 'REV']))
      an <- an[i,, drop=FALSE]
      rownames(an) <- sub('  (Factor+Higher Order Factors)', '',
                          rownames(an), fixed=TRUE)

      if(isbase) {
        xlim <- range(an[, 1:3])
        dotchart2(an[, 'REV'], xlab=xlab, pch=pch, xlim=xlim, ...)
        dotchart2(an[, 'Lower'], pch=91, add=TRUE)
        dotchart2(an[, 'Upper'], pch=93, add=TRUE)
        return(invisible(an))
      }
      p <- dotchartpl(an[, 'REV'], major=rownames(an),
                      lower=an[,'Lower'], upper=an[,'Upper'],
                      xlab=xlab,
                      limitstracename='HPD Interval',
                      width=width, height=height)
      return(p)
      }

    if(what %in% c("partial R2", "remaining R2", "proportion R2")) {
      if("Partial SS" %nin% colnames(x))
        stop('to plot R2 you must have an ols model and must not have specified ss=FALSE to anova')

      sse <- x ['ERROR', 'Partial SS']
      ssr <- x ['TOTAL', 'Partial SS']
      pss <- an[, 'Partial SS']
      sst <- sse + ssr
    }

    dof <- an[, 'd.f.']
    P   <- an[, 'P']

    if(any(colnames(an) == 'F')) {
      chisq    <- an[, 'F'] * dof
      totchisq <- x['TOTAL', 'F'] * x['TOTAL', 'd.f.']
    }
    else {
      chisq    <- an[, 'Chi-Square']
      totchisq <- x['TOTAL', 'Chi-Square']
    }

    w <- switch(what,
                chisq           = chisq,
                chisqminusdf    = chisq - dof,
                aic             = chisq - 2 * dof,
                P               = P,
                "partial R2"    = pss / sst,
                "remaining R2"  = (ssr - pss) / sst,
                "proportion R2" = pss / ssr,
                "proportion chisq" = chisq / totchisq)

    if(missing(newnames))
      newnames <- sedit(names(w),"  (Factor+Higher Order Factors)", "")

    names(w) <- newnames
    is <- switch(sort,
                descending =  order(-w),
                ascending  =  order( w),
                none       =  1 : length(w))
    w     <- w [is]
    an    <- an[is,, drop=FALSE ]
    chisq <- chisq[is]
    dof   <- dof[is]
    P     <- P[is]

    if(pl) {
      auxtitle <- auxdata <- NULL
      fn <- function(x, right) {
        m <- max(abs(x), na.rm=TRUE)
        left <- max(floor(log10(m)) + 1, 1)
        nFm(x, left, right)
      }

      if(any(c('partial R2', 'remaining R2') %in% margin)) {
        if("Partial SS" %nin% colnames(x))
          stop('to show R2 you must have an ols model and must not have specified ss=FALSE to anova')
        sse <- x['ERROR', 'Partial SS']
        ssr <- x['TOTAL', 'Partial SS']
        sst <- sse + ssr
        pss <- an[, 'Partial SS']
      }

      if(length(margin))
        for(marg in margin) {
          aux <-
            if(isbase)
              switch(marg,
                     chisq = list('chi^2', fn(chisq, 1)),
                     'proportion chisq' =
                       list('Proportion~chi^2', fn(chisq / totchisq, 2)),
                     'd.f.' = list('d.f.', fn(dof, 0)),
                     P = list('P', fn(P, 4)),
                     'partial R2' = list('Partial~R^2',       fn(pss / sst, 2)),
                     'proportion R2' = list('Proportion~R^2', fn(pss / ssr, 2)))
            else
              switch(marg,
                     chisq = paste(htmlSpecs$chisq(dof), fn(chisq, 1)),
                     'proportion chisq' =
                       paste0('Proportion ', schisq, '=',
                             fn(chisq / totchisq, 2)),
                     'd.f.' = paste('d.f.=', fn(dof, 0)),
                     P = paste('P=', fn(P, 4)),
                     'partial R2' = paste('Partial R<sup>2</sup>=',
                                          fn(pss / sst, 2)),
                     'proportion R2' = paste('Proportion R<sup>2</sup>=',
                                             fn(pss / ssr, 2)))

          if(isbase) {
            if(length(auxtitle))
              auxtitle <- paste(auxtitle, aux[[1]], sep='~~')
            else auxtitle <- aux[[1]]
            if(length(auxdata))
              auxdata  <- paste(auxdata,  aux[[2]], sep='  ')
            else auxdata  <- aux[[2]]
          } else
            auxdata <- if(length(auxdata))
                         paste(auxdata, aux, sep=paste0(nbsp,nbsp))
                       else
                         aux
      }
      ## convert to expression if not using plotly
      if(length(auxtitle) && isbase) auxtitle <- parse(text = auxtitle)

      dc <- if(isbase) dotchart3 else dotchartp
      if(length(trans)) {
        nan <- names(w)
        w <- pmax(0, w)
        pan <- pretty(w, n=ntrans)
        tan <- trans(w); names(tan) <- nan
        p <- dc(tan, xlab=xlab, pch=pch,
                axisat=trans(pan), axislabels=pan,
                auxtitle=auxtitle, auxdata=auxdata, auxwhere='hover',
                height=height, width=width, ...)
      } else p <- dc(w, xlab=xlab, pch=pch,
                     auxtitle=auxtitle, auxdata=auxdata, auxwhere='hover',
                     height=height, width=width, ...)
    }
    if(isbase) invisible(w) else p
  }

## -----------------------------------------------------------------------
## Detailed notes on extending latex.anova.rms() to typst, referenced by
## the short [N] tags inline above. Most of this is simply building out
## Typst coverage for the first time -- part of this project's goal from
## the start, not corrections to anything broken.
##
## [1] rowl/resp translation: previously `if(! html) rowl <-
##     latexTranslate(rowl)` applied LaTeX-specific escaping to any
##     non-html language, since typst wasn't part of the original
##     design. Now lang-aware: latexTranslate() for latex,
##     typstTranslate() for typst, untouched for html/plain. Without
##     this, LaTeX escape sequences (e.g. "\%", "\_") would appear
##     literally in Typst output rather than being escaped correctly for
##     Typst's own syntax. Same extension applied to resp (the response
##     variable name used in the caption).
##
## [2] Interaction symbol (*) substitution: extending this to typst
##     needs the search pattern to match what typstTranslate() has
##     already produced by that point, not a bare '*'. typstTranslate()
##     escapes bare * to \* unconditionally (needed in general, since *
##     is Typst's bold/emphasis shorthand) -- so with a bare '*' pattern,
##     fixed=TRUE would match the * inside \*, leaving the escaping
##     backslash orphaned and producing a broken token like "\$times$".
##     Same construct already handled this way in prModFit's coefficient
##     table row-name handling, kept consistent here. latex/html's
##     translators never touch *, so they still use the bare pattern.
##
## [3] Table-making, rms_tiny_table(): replaces the previous
##     htmlTable::htmlTable()/Hmisc::latex() split, which had no typst
##     handling at all. Structured as one format-agnostic core builder
##     (rms_tiny_table_core) plus a separate, named finalize and caption
##     function per language (rms_tiny_table_finalize_latex/html/typst,
##     rms_tiny_table_caption_latex/html/typst), so each format's
##     handling is self-contained and can be extended independently
##     without touching the others -- this is also intended to be
##     reused by summary.rms and validate.ols, which produce
##     structurally similar tables.
##
##     table.env's former role (choosing between a numbered/captioned
##     LaTeX float vs. a plain caption) is no longer used --
##     environment_table=FALSE is applied unconditionally, matching the
##     approach used throughout this project to avoid LaTeX's
##     auto-numbered "Table N:" captioning; the caption is always
##     emitted as a plain bold, centered heading instead via
##     rms_tiny_table_caption_latex/html/typst.
##
##     specs$lspace (used in the leading-blank-preservation logic,
##     unchanged from the original) is confirmed defined in
##     markupSpecs$typst as "#h(0.5em)".
##
## [4] Chi-square column name: `sn[sn=='Chi-Square'] <- specs$chisq(add='')`
##     followed by `ifelse(sn %nin% c(...), math(sn), sn)` applied math()
##     a second time to the same entry -- by the time the second line's
##     %nin% check ran, sn no longer held "Chi-Square" for that column
##     (it had already been overwritten with the chisq()-formatted
##     value), so the check couldn't exclude it and math() ran again.
##     Under typst this double-wrapped already-$-wrapped content,
##     surfacing as literal unrendered text. See note [5]: the original
##     double-application turned out to be load-bearing for latex,
##     where it was the only thing providing $ wrapping at all -- fixed
##     properly there, not just removed.
##
## [5] specs$chisq()'s $-wrapping turns out to be format-inconsistent:
##     typst's version already returns fully-wrapped math ("$chi^2$"),
##     while latex's returns bare, unwrapped "\chi^{2}", relying on
##     something downstream to add the $ delimiters -- which, before
##     note [4]'s fix, was exactly the "double" math() application being
##     removed there. Removing it without accounting for this asymmetry
##     fixed typst but broke latex outright: bare \chi (a math-mode-only
##     command) outside $ $ produces "Missing $ inserted" at LaTeX
##     compile time, confirmed directly from a real compile error.
##     Fixed by checking the actual returned value rather than assuming
##     either format's behavior: if specs$chisq()'s output already
##     starts with "$", use it as-is; otherwise wrap it with math().
##     This is format-agnostic -- it adapts to whatever each language's
##     chisq() actually returns rather than hard-coding an assumption
##     per format, so it should hold even if a given format's wrapping
##     convention isn't fully known or changes later.
##
## [6] rms_tiny_table_core/rms_tiny_table extended with several new,
##     backward-compatible optional parameters (defaults preserve prior
##     behavior exactly, so anova.rms's and summary.rms's existing calls
##     are unaffected): show_rownames=FALSE omits the row-label column
##     entirely, needed by validate.ols for its marker/frequency tables,
##     which have no meaningful row names (the original code used
##     rowname=NULL/rnames=FALSE for these). format_tt(linebreak="\n")
##     is now applied unconditionally to headers, enabling multi-line
##     column headers (e.g. "Original\nSample") -- harmless for headers
##     without embedded newlines, and reusing the exact mechanism
##     already confirmed working for prStats' multi-line headers.
##
## [7] Refactoring pass across prModFit/prStats (rmsMisc.s), prompted by
##     a review specifically looking for redundancy across the model-fit
##     print/latex methods and prModFit, once all of anova.rms/
##     summary.rms/validate.ols/prModFit/prStats existed side by side:
##
##     coef_tiny_table (rmsMisc.s, prModFit's coefficient table) was
##     deleted entirely -- comparing it directly against
##     rms_tiny_table_core, it was a strict subset (same data.frame/
##     escape/alignment/centering/striping construction, built before
##     rms_tiny_table existed and never reconciled with it afterward).
##     Its one genuinely new contribution, the discovery that cph/Gls/
##     psm's real output wraps prModFit's content in an outer #block[ of
##     unknown origin with no definite width (breaking simple
##     #align(center) alone), was already independently incorporated
##     into rms_tiny_table_finalize_typst's #block(width:100%) wrap, so
##     nothing was lost by removing it. prModFit's one call site now
##     calls rms_tiny_table(U, lang) directly.
##
##     stats_tiny_table (rmsMisc.s, prStats' final table) was NOT a
##     simple subset -- it has real, hard-won features rms_tiny_table
##     didn't: no row-label column at all (every column holds
##     dual-justified "label #h(1fr) value" content via hfill, not a
##     separate label + right-aligned value column), a bold header with
##     a bottom rule, vertical column borders, column-count-scaled
##     width, and format-specific spacing tuning (LaTeX colsep/rowsep,
##     Typst fontsize reduction) -- each added after a real bug surfaced
##     during testing. Extended rms_tiny_table's parameters (align_data,
##     bold_header, header_rule, vertical_borders, width, latex_inner --
##     see note [6] immediately above) to cover all of these, all
##     defaulting to prior behavior exactly. stats_tiny_table itself is
##     now a thin wrapper: applies its own no-hyphenation header wrapping
##     (genuinely specific to this table, kept local) and computes its
##     format-specific width/fontsize/inner values, then makes one call
##     to rms_tiny_table() with those as overrides.
##
## [8] Real regression, caught immediately via a real LaTeX compile
##     error ("Missing $ inserted") on the very next test after note
##     [7]'s refactor: rms_tiny_table_core's show_rownames=FALSE branch
##     added an unclass(d)/as.data.frame(..., check.names=FALSE)
##     round-trip that the original stats_tiny_table() never had --
##     it called tt() directly on its already-built data.frame, no
##     conversion step at all. validate.ols's varin/tkept tables also
##     go through this same show_rownames=FALSE branch and were already
##     confirmed working, so the round-trip itself isn't universally
##     broken; the difference is that stats_tiny_table's column headers
##     are multi-line and carry embedded Typst/LaTeX markup (built via
##     colnames(d) <- labels), unlike validate.ols's simple column
##     names, and something about that combination interacting with the
##     round-trip corrupted the output -- confirmed directly from the
##     generated .tex file, which showed an entire column's data
##     collapsed into one cell as literal, deparsed R vector syntax
##     (c("Obs~\\hfill 960", ...)) rather than being used as actual
##     per-row values. Rather than fully diagnose the exact
##     unclass()/as.data.frame() interaction, fixed by skipping the
##     round-trip entirely when d is already a data.frame -- true for
##     every current caller, so this sidesteps the issue without
##     depending on understanding it completely.
##
## [9] align_data extended to accept a per-column vector, not just a
##     single value applied uniformly (the prior, still-default
##     behavior). Prompted by a real crash in processMI.r's prmiInfo():
##     it had been converting a "Test" label column into data.frame
##     rownames (rownames(m) <- m$Test) to get that column
##     left-aligned via show_rownames=TRUE, but data.frame rownames must
##     be unique, and prmiInfo's Test labels legitimately are not --
##     multiple different predictors can each have their own
##     "Nonlinear"/"All Interactions" sub-test row, sharing the same
##     display text for the test type while representing genuinely
##     different predictors. This was a real, confirmed regression: the
##     original prmiInfo (before being unified with this shared helper)
##     called htmlTable::htmlTable(m, ..., rnames=FALSE) specifically to
##     avoid this exact constraint, keeping Test as an ordinary data
##     column rather than rownames. Fixed at the source in prmiInfo
##     (kept as an ordinary column, show_rownames=FALSE now), and this
##     align_data extension lets it still get the same per-column
##     layout as the original (Test left-aligned, the three numeric
##     columns right-aligned) without needing rownames at all. A vector
##     align_data loops style_tt() once per column rather than trusting
##     element-wise recycling of align across a vector j (unconfirmed
##     for tinytable's style_tt()) -- the same safe, already-proven
##     per-column style_tt() pattern used elsewhere in this function.
## -----------------------------------------------------------------------
