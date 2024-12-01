#main.effect=F to suppress printing main effects when the factor in
#question is involved in any interaction.

anova.rms <- function(object, ..., main.effect=FALSE, tol=1e-13, 
                      test=c('F', 'Chisq', 'LR'),
                      india=TRUE, indnl=TRUE, ss=TRUE,
                      vnames=c('names', 'labels'),
                      posterior.summary=c('mean', 'median', 'mode'),
                      ns=500, cint=0.95, fitargs=NULL) {

  misstest <- missing(test)
  test     <- match.arg(test)

  ava <- if(test == 'LR') function(idx) LRchunktest(object, idx, tol=tol, fitargs=fitargs)
            else function(idx) {
              chisq <- coef[idx] %*% solvet(cov[idx, idx], coef[idx], tol=tol)
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
    
    if(! html) rowl <- latexTranslate(rowl)

    specs <- markupSpecs[[lang]]
    bold  <- specs$bold
    math  <- specs$math
    
    ## Translate interaction symbol (*) to times symbol
    ## rowl <- gsub('\\*', specs$times, rowl)   # changed * to $times$
    rowl <- gsub('*', specs$times, rowl, fixed=TRUE)
  
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
    sn[sn=='Chi-Square'] <- specs$chisq(add='')
    names(dstats) <- ifelse(sn %nin% c('d.f.','MS','Partial SS'),
                            math(sn), sn)

    resp <- as.character(attr(object, 'formula')[2])
    if(! html) resp <- latexTranslate(resp)

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

    if(html) {
      al <- rep('r', length(sn))
      fshead <- rep(paste0('font-size:', fontsize, 'em;'), ncol(dstats))
      fscell <- rep('padding-left:2ex;',                   ncol(dstats))
      w <- htmlTable::htmlTable(dstats, caption=caption,
                                css.table=fshead,
                                css.cell =fscell,
                                align=al, align.header=al,
                                rowlabel='', escape.html=FALSE)
      rendHTML(w)
    }
    else {
        latex(dstats, title=title,
              caption    = if(table.env) caption else NULL,
              insert.top = if(length(caption) && ! table.env)
                             paste0('\\Needspace{2in}\n', caption),
              rowlabel="", col.just=rep('r',length(sn)), table.env=table.env, ...)
      }
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
