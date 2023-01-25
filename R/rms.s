# Design  FEH 1Aug90, re-written 21Oct91
#
# Augments S formula language to include:
# 
#	name	- name[i] = name of ith original variable in x
#	label	- label[i] = label of ith original variable (=name if none)
#	assume	- assume(original x)
#	assume.code - coded version of assume (1-11, 9=added interaction)
#	parms	- parms(original x)
#		  for interaction effects parms[[i]] is a matrix with dim
#		  3 x (1+# interaction terms).  First element in pair
#		  is 1 if first factor is represented as an expanded
#		  non-linear term, 0 otherwise (this applies to polynomial,
#		  lspline, rcspline, scored).  Second element applies to
#		  second factor in interaction effect.  Third element
#		  applies to third factor (0 if second order interaction)
#		  First column contains factor numbers involved in interaction.
#	limits  - limits(original x)
#	values	- For continuous variables with <=10 unique values, is
#		  vector of values.  NULL otherwise.
#	interactions - 3 x k matrix of factor numbers
#
# Cannot have interactions between two stratification factors. 
#
#

  
Design <- function(mf, formula=NULL, specials=NULL,
                   allow.offset=TRUE, intercept=1) {

  debug <- length(.Options$rmsdebug) && .Options$rmsdebug
  if(debug) {
    cat('--------------------------------------------------\n')
    prn(list(names(mf), formula, specials))
    }

  if(length(formula)) {
    Terms <- terms(formula, specials=specials, data=mf)
    attr(mf, 'terms') <- Terms
    }
  else Terms <- attr(mf, 'terms')

  Terms.orig  <- Terms
  Term.labels <- attr(Terms, 'term.labels')
  ## offsets are not included anywhere in terms even though they are
  ## in the model frame
  response.pres <- attr(Terms, 'response') > 0

  ## Function to construct the colname names that model.matrix will
  ## create.  This is primarily used to subset the columns of the model
  ## matrix to get rid of terms involving strat main effects and to get
  ## rid of interaction terms involving non-reference values

  mmnames <- function(assume.code, rmstrans.names, term.label, iaspecial,
                      class) {
    if(debug) {
      prn(assume.code); prn(rmstrans.names); prn(term.label); prn(iaspecial); prn(class) }
    ## Don't let >=i be translated to >i:
    rmst <- gsub('>=', '>>', rmstrans.names)
    ## Don't let <=i be translated to <i:
    rmst <- gsub('<=', '<<', rmst)
    ## Don't let == be translated to blank
    rmst <- gsub('==', '@EQ@', rmst)
    if(debug) {prn(class);prn(assume.code);prn(term.label)}
    w <- if(assume.code == 1)
           ifelse(class == 'logical', paste0(term.label, 'TRUE'),
                  term.label)
    else if(length(iaspecial) && iaspecial) term.label
    else if(assume.code == 5) {
      ## Handle explicit catg() as well as implicit
      ## Use sub to only change the first occurrence; value
      ## labels may contain = also
      if(substring(term.label, 1, 5) == 'catg(')
        paste0(term.label, sub('.*=', '', rmst))
      else sub('=', '', rmst)
      }
    else if(assume.code == 8)
      paste0(term.label, sub('.*=', '', rmst))
    else if(assume.code >= 10)   # was == 10
      if(length(rmst) > 1) gsub('\\[', '', gsub('\\]', '', rmst)) else
             term.label
    else paste0(term.label, rmst)
    w <- gsub('>>',   '>=', w)
    w <- gsub('<<',   '<=', w)
    w <- gsub('@EQ@', '==', w)
    alt <- if(assume.code >= 10)   # was == 10
             if(length(rmst) > 1)
               paste0(term.label, rmstrans.names) else term.label
    else w
    ## Alternate names to try - handles case where model is fitted on a
    ## previous fit$x matrix
    attr(w, 'alt') <- alt
    w
  }

  offs <- model.offset(mf)
  iscluster <- if(length(Term.labels))
                 substring(Term.labels, 1, 8) == 'cluster('  else FALSE
  istime <- if(length(Term.labels))
              substring(Term.labels, 1, 6) == 'aTime('  else FALSE

  ## Handle cluster() and aTime()
  ## Save right hand side of formula less cluster() and time() terms
  sformula <- formula(Terms)
  if(any(iscluster)) sformula <- removeFormulaTerms(sformula, 'cluster')
  if(any(istime))    sformula <- removeFormulaTerms(sformula, 'aTime')

  if(any(iscluster)) {
    clustername <- Term.labels[iscluster]
    cluster     <- mf[[clustername]]
    mf[[clustername]] <- NULL
    Terms       <- Terms[! iscluster]
    Term.labels <- Term.labels[! iscluster]
    if(any(istime)) istime <- if(length(Term.labels))
     substring(Term.labels, 1, 6) == 'aTime('  else FALSE
  }
  else {cluster <- clustername <- NULL}

  if(any(istime)) {
    timename <- Term.labels[istime]
    time     <- mf[[timename]]
    mf[[timename]] <- NULL
    Terms       <- Terms[! istime]
    Term.labels <- Term.labels[! istime]
  }
  else {time <- timename <- NULL}
  

  ioffset <- integer(0)
  if(length(offs)) {
    if(! allow.offset) stop("offset variable not allowed in formula")
    ## first case below is with offset= in fit call, 2nd with offset(var)
    ioffset <- which(names(mf) == '(offset)' |
                       substring(names(mf), 1, 7) == 'offset(')
    if(! any(ioffset)) stop('program logic error 1')
  }

  ## For some reason, model frame sometimes has a blank name if using %ia%

  namx <- names(mf)
  if(any(namx == "")) {
    namx <- names(mf) <- c(namx[1], Term.labels)
    dimnames(mf)[[2]] <- namx
    dimnames(attr(Terms, "factors"))[[1]] <- namx
  }

  wts <- if(any(namx == '(weights)'))(1 : length(namx))[namx == '(weights)']
         else 0

  if(debug) prn(llist(names(mf), ioffset, response.pres, wts))

  coluse <- setdiff(1 : ncol(mf), c(ioffset, 1 * response.pres, wts))

  inner.name <- if(length(Terms) > 0) unique(var.inner(Terms))
  ## Handles case where a function has two arguments that are names,
  ## e.g. rcs(x,knots) -> want x only
  ## Note: these exclude interaction terms and %ia% terms
  
  factors <- attr(Terms, "factors")
  if(length(factors) && response.pres) factors <- factors[-1, , drop=FALSE]

  attr(Terms, "intercept") <- intercept

  ## name:   nice version of design matrix column names constructed here
  ## mmname: column names that model.matrix will later create
  ##         Used for later keeping only those columns that don't pertain
  ##         to strat main effects or to strat interactions involving
  ##         non-reference categories
  
  fname <- flabel <- name <- mmname <- strt <- asm <- len <- 
    fname.incl.dup <- ia <- funits <- NULL
  parm <- nonlinear <- tex <- limits <- values <- list()

  scol      <- 1
  colnam    <- mmcolnam <- mmcolnamalt <- list()
  Altcolnam <- NULL

  XDATADIST <- .Options$datadist
  if(length(XDATADIST)) {
    if(inherits(XDATADIST, 'datadist')) datadist <- XDATADIST
       else {
         if(! exists(XDATADIST))
           stop(paste("dataset", XDATADIST,
                      "not found for options(datadist=)"))
         datadist <- eval(as.name(XDATADIST))
         }
    Limits   <- datadist$limits
    Limnames <- dimnames(Limits)[[2]]
  }
  
  nc <- 0

  options(Design.attr=NULL, TEMPORARY=FALSE)
  ##Used internally by asis, rcs, ...

  anyfactors <- length(coluse) > 0
  i1.noia <- 0
  if(length(Term.labels) < length(coluse))
    stop(paste('program logic error tl\nTerm.labels:',
               paste(Term.labels, collapse=', '), '\ncoluse:',
               paste(coluse, collapse=', ')))
  it <- 0
  if(anyfactors) for(i in coluse) {
    if(i  != wts) {
      i1 <- i - response.pres
      xi <- mf[[i]]
      cls <- rev(class(xi))[1]
      z <- attributes(xi)
      assu <- z$assume.code
      if(! length(assu) || assu != 9) i1.noia <- i1.noia + 1
      if(! length(assu)) {
        ## Not processed w/asis,et
        nam <- inner.name[i1.noia]
        lab <- attr(xi, "label")
        ord <- is.ordered(xi) && all.is.numeric(levels(xi))
        if(! length(lab) || lab == "") lab <- nam
        if(ord) {
          xi <- scored(xi, name=nam, label=lab)
          attr(mf[, i], "contrasts") <- attr(xi, "contrasts")
        }
        else if(is.character(xi) | is.factor(xi)) {
          if(is.ordered(xi) &&
             .Options$contrasts[2] != 'contr.treatment')
            stop(paste('Variable', nam,
                       'is an ordered factor with non-numeric levels.\n',
                       'You should set options(contrasts=c("contr.treatment", "contr.treatment"))\nor rms will not work properly.'))
          xi <- catg(xi, name=nam, label=lab)
        }
        else if(is.matrix(xi)) xi <- matrx(xi, name=nam, label=lab)
        else xi <- asis(xi, name=nam, label=lab)
        z <- c(z, attributes(xi))
      }

      za    <- z$assume.code
      zname <- z$name

      fname.incl.dup <- c(fname.incl.dup, zname)
      if(! length(fname) || ! any(fname == zname)) { # unique factor
        nc <- nc + 1
        fname  <- c(fname,  zname)
        flabel <- c(flabel, z$label)
        asm <- c(asm, za)
        colnam[[i1]] <- z$colnames
        it <- it + 1
        mmn <- mmnames(za, colnam[[i1]], Term.labels[it], z$iaspecial, cls)
        mmcolnam[[i1]] <- mmn
        alt <- attr(mmn, 'alt')
        mmcolnamalt[[i1]] <- alt
        if(debug) prn(c(mmn, alt))
        if(za != 8 && length(colnam)) {
          name   <- c(name,   colnam  [[i1]])
          mmname <- c(mmname, mmcolnam[[i1]])
          Altcolnam <- c(Altcolnam, alt)
        }
        if(za != 9) {
          funits <- c(funits, if(length(z$units))z$units else '')
          if(length(z$parms)) parm[[zname]] <- z$parms
          if(length(XDATADIST)) {
            limits[[zname]] <- if(any(Limnames == zname)) {
              j <- match(zname, Limnames, 0) #require EXACT match
              Limits[, j[j > 0]]
            }
            else rep(NA, 7)
            j <- match(zname, names(datadist$values), 0)
            if(j > 0) {
              values[[zname]] <- datadist$values[[j]]
              l1 <- levels(xi); l2 <- datadist$values[[j]]
              if(length(l1) && ((length(l1)  != length(l2)) ||
                                any(sort(l1) != sort(l2))))
                warning(paste('Variable', zname, 'has levels', paste(l1, collapse=' '),
                              'which do not match levels given to datadist (',
                              paste(l2, collapse=' '), '). datadist values ignored.'))
              values[[zname]] <- l1
            }
          }
        }
        
        if(length(nonl <- z$nonlinear)) nonlinear[[zname]] <- nonl

        if(length(tx <- z$tex))         tex[[zname]]       <- tx

        if(za == 9) {
          iia <- match(z$ia, fname)
          if(any(is.na(iia)))stop(paste(paste(z$ia, collapse=" "),
                "cannot be used in %ia% since not listed as main effect"))
          ia <- cbind(ia, c(iia, 0))
          parms <- rbind(z$parms, 0)
          parms[, 1] <- c(iia, 0)
          if(length(parms)) parm[[zname]] <- parms
        }
      }
      nrows <- if(is.matrix(xi)) nrow(xi) else length(xi)
    }
  }
  
  ##Save list of which factors where %ia% interactions
  ## (before adding automatic ias)

  which.ia <- (1 : length(asm))[asm == 9]

  ##Add automatically created interaction terms
  
  if(anyfactors) {
    nrf <- if(! length(factors)) 0 else nrow(factors)

    if(length(factors)) for(i in 1 : ncol(factors)) {
      f <- factors[, i]
      j <- (1 : length(f))[f > 0]
      nia <- length(j)
      if(nia > 1) {
        fn <- fname.incl.dup[j]
        jf <- match(fn, fname.incl.dup)
        if(any(is.na(jf))) stop("program logic error 2")
        nc <- nc + 1
        asm <- c(asm, 9)
        if(nia == 2) ialab <- paste(fn[1], "*", fn[2])
        else if(nia == 3)ialab <- paste(fn[1], "*", fn[2], "*", fn[3])
        else stop("interaction term not second or third order")
        fname  <- c(fname,  ialab)
        flabel <- c(flabel, ialab)
        if(sum(asm[jf] == 8) > 1)
          stop("cannot have interaction between two strata factors")
        nn <- mmnn <- mmnnalt <- list()
        for(k in 1 : nia) {
          if(asm[jf[k]] == 5 | asm[jf[k]] == 8)
            nn[[k]] <- paste0(fn[k], "=", parm[[fname[jf[k]]]][-1])
          else if(asm[jf[k]] == 7) {
            nn[[k]] <- c(fn[k],
                         paste0(fn[k], "=",
                               parm[[fname[jf[k]]]][c(-1, -2)]))
          }
          else nn[[k]] <- colnam[[jf[k]]]
          mmnn[[k]]    <- mmcolnam[[jf[k]]]
          mmnnalt[[k]] <- mmcolnamalt[[jf[k]]]
        }
        if(nia == 2) {nn[[3]] <- mmnn[[3]] <- mmnnalt[[3]] <- ""}
        parms <- jf
        if(length(jf) == 2) parms <- c(parms, 0)
        nonlin <- NULL
        nl1 <- nonlinear[[fname[jf[1]]]]
        nl2 <- nonlinear[[fname[jf[2]]]]
        ## Strata factors don't have nonlinear duplicated for # levels - 1
        if(asm[jf[1]] == 8)
          nl1 <- rep(FALSE, length(parm[[fname[jf[1]]]]) - 1)
        if(asm[jf[2]] == 8)
          nl2 <- rep(FALSE, length(parm[[fname[jf[2]]]]) - 1)
        if(nia == 2) nl3 <- FALSE
        else if(asm[jf[3]] == 8)
          nl3 <- rep(FALSE, length(parm[[fname[jf[3]]]]) - 1)
        else nl3 <- nonlinear[[fname[jf[3]]]]
        n1 <- nn[[1]]
        n2 <- nn[[2]]
        n3 <- nn[[3]]
        mmn1 <- mmnn[[1]]
        mmn2 <- mmnn[[2]]
        mmn3 <- mmnn[[3]]
        mmnalt1 <- mmnnalt[[1]]
        mmnalt2 <- mmnnalt[[2]]
        mmnalt3 <- mmnnalt[[3]]
        
        ## model.matrix makes auto-products move first variable fastest, etc.
        for(j3 in 1 : length(n3)) {
          for(j2 in 1 : length(n2)) {
            for(j1 in 1 : length(n1)) {
              parms <- cbind(parms, c(nl1[j1], nl2[j2], nl3[j3]))
              nonlin <- c(nonlin, nl1[j1] | nl2[j2] | nl3[j3])
              name <- c(name,
                        if(nia == 2) paste(n1[j1], "*", n2[j2])
                        else paste(n1[j1], "*", n2[j2], "*", n3[j3]))
              mmname <- c(mmname,
                          if(nia == 2) paste0(mmn1[j1], ':', mmn2[j2])
                          else paste0(mmn1[j1], ':', mmn2[j2], ':', mmn3[j3]))
              Altcolnam <- c(Altcolnam,
                             if(nia == 2) paste0(mmnalt1[j1], ':', mmnalt2[j2])
                             else paste0(mmnalt1[j1], ':', mmnalt2[j2], ':',
                                        mmnalt3[j3]))
            }
          }
        }
        
        ## If was 2-way interaction and one of the factors was restricted %ia%,
        ## adjust indicators
        k <- match(jf, which.ia, 0)
        if(any(k > 0)) {
          if(nia == 3)
            stop("cannot have 2-way interaction with an %ia% interaction")
          k <- jf[k > 0]
          wparm <- parms[, 1] == k; wparm[3] <- TRUE
          parms[wparm,] <- parm[[fname[k]]][1 : 2,, drop=FALSE]
          jf <- parms[, 1]
          nonlin <- apply(parms, 2, any)[-1]
        }
        
        if(length(jf) == 2) {jf <- c(jf, 0); parms[3, ] <- 0}
        ia <- cbind(ia, jf)
        if(length(parms)) parm[[ialab]] <- parms
        if(length(nonlin)) nonlinear[[ialab]] <- nonlin
        
      }
    }
  }
  
  if(anyfactors) {
    if(length(XDATADIST))
      limits <- structure(limits,
                          row.names=c("Low:effect", "Adjust to",
                            "High:effect",
                            "Low:prediction", "High:prediction",
                            "Low", "High"),
                          class="data.frame")
    ##data.frame converts variables always NA to factor!
    
    if(length(funits) != sum(asm != 9)) warning('program logic warning 1')
    else names(funits) <- fname[asm != 9]

    attr(mmname, 'alt') <- if(! all(Altcolnam == mmname)) Altcolnam
    if(any(duplicated(mmname)))
      stop(paste0('duplicated column name in design matrix:',
                 paste(mmname[duplicated(mmname)], collapse=' '),
                 '\nMost likely caused by a variable name concatenated to a factor level\nbeing the same is the name of another variable.'))
    atr <- list(name=fname, label=flabel, units=funits,
                colnames=name, mmcolnames=mmname,
                assume=c("asis", "polynomial", "lspline", "rcspline",
                  "category", "","scored", "strata", "interaction",
                  "matrix", "gTrans")[asm],
                assume.code=as.integer(asm), parms=parm, limits=limits,
                values=values, nonlinear=nonlinear, tex=tex,
                interactions=if(length(ia)) structure(ia, dimnames=NULL))
    
    nact <- attr(mf, 'na.action')
    if(length(nact) && length(nmiss <- nact$nmiss)) {
      jia <- grep('%ia%', names(nmiss),  fixed=TRUE)
      if(length(jia)) nmiss <- nmiss[-jia]
      jz <- which(names(nmiss) != '(weights)' &
                    ! grepl('offset\\(',  names(nmiss)) &
                    names(nmiss) != '(offset)' &
                    ! grepl('cluster\\(', names(nmiss)) &
                    ! grepl('aTime\\(',   names(nmiss)))
      if(response.pres) jz <- jz[jz > 1]
      names(nmiss)[jz] <- fname[asm != 9]
      attr(mf, 'na.action')$nmiss <- nmiss
    }
  }
  
  else atr <- list(name=NULL, assume=NULL, assume.code=NULL, parms=NULL)
  
  attr(mf, 'Design')     <- atr
  attr(mf, 'terms')      <- Terms
  attr(mf, 'sformula') <- sformula
  if(length(cluster)) {
    attr(mf, 'cluster') <- cluster
    attr(mf, 'clustername') <- var.inner(as.formula(paste0('~', clustername)))
  }
  if(length(time)) {
    attr(mf, 'time') <- time
    attr(mf, 'timename') <- var.inner(as.formula(paste0('~', timename)))
  }
  
  if(length(offs))    attr(mf, 'offset')  <- offs
  mf
}

modelData <- function(data=environment(formula), formula, formula2=NULL,
                      weights=NULL, subset=NULL, na.action=na.delete,
                      dotexpand=TRUE, callenv=parent.frame(n=2)) {

  ## calibrate.cph etc. uses a matrix, even if only one column
  ismat <- function(z) {
    cl <- class(z)
    ('matrix' %in% cl) && ('rms' %nin% cl) ## && ncol(z) > 1
    }
  
  ## Get a list of all variables in either formula
  ## This is for innermost variables, e.g. Surv(a,b) will produce a,b
  v1 <- all.vars(formula)
  v2 <- all.vars(formula2)
  V  <- unique(c(v1, v2))

  edata <- is.environment(data)

  rhsdot <- length(v1) == 2 && v1[2] == '.'
  if(rhsdot && edata)
    stop('may not specify ~ . in formula when data= is absent')

  if(edata) {
    env  <- data
    data <- list()
    for(v in V) {
      xv <- env[[v]]
      if(is.factor(xv)) xv <- xv[, drop=TRUE]
      ## Note: Surv() has class 'Surv' without class 'matrix'
      ## This keeps columns together by calling as.data.frame.rms
      if(ismat(xv)) class(xv) <- unique(c('rms', class(xv)))
      data[[v]] <- xv
    }
    ## Any variables whose length is not equal to the maximum length over
    ## all variables mentioned in the formulas remain in the original
    ## environment and will be found in the later eval()
    ## E.g. rcs(x, knots) where knots is a separate variable
    n <- sapply(data, NROW)
    if(! length(n)) stop('no data found')
    if(diff(range(n)) != 0) data <- data[which(n == max(n))]
    ## Watch out: if a variable in data has dimnames[[2]], as.data.frame
    ## uses that as the new variable name even if the variable already
    ## had a name in the list.  This is why a 1-column matrix is kept
    ## as a matrix in the ismat function above
    data <- as.data.frame(data)
  }   # end if(edata)
  ## Can't do else data[V] here as formula may have e.g. Surv(time,event)
  ## and hasn't been evaluated yet, where data has time and event
  if(length(weights)) data$`(weights)` <- weights

  if(length(subset)) data <- data[subset, ]

  ## Make sure that the second formula doesn't create any NAs on
  ## observations that didn't already have an NA for variables in main formula
  if(length(formula2)) {
    i <- ! complete.cases(data[intersect(names(data), v1)])
    j <- ! complete.cases(data[intersect(names(data), v2)])
    ## Check: longer object length not mult of shorter:  ???
    if(any(j & ! i))
      stop('A variable in the second formula was missing on an observation that was not missing on any variable in main formula')
  }

  noexpand <- rhsdot & ! dotexpand
  deparse2 <- function(x)   # from stats
    paste(deparse(x, width.cutoff = 500L, backtick = !is.symbol(x) && 
                  is.language(x)), collapse = " ")

  processdata <- function(formula, data) {
    if(noexpand) {   # no RHS variables to be used
      predvars <- formula[[2]]
      varnames <- deparse(predvars)
      if(length(weights)) {
        predvars[[2]] <- as.name('(weights)')
        varnames      <- c(varnames, '(weights)')
      }
    }
    else {
      Terms    <- terms(formula, data=data, specials=NULL)
      vars     <- attr(Terms, 'variables')
      predvars <- attr(Terms, 'predvars')
      if( ! length(predvars)) predvars <- vars
      if(length(weights))
        predvars[[length(predvars) + 1]] <- as.name('(weights)')
    }
    varnames <- vapply(predvars, deparse2, " ")[-1L]
    rnames <- rownames(data)     # 2023-01-25
    data <- if(edata) eval(predvars, data, env)
            else
              eval(predvars, data, callenv)

    if(is.matrix(data)) data <- data.frame(data)  # e.g. Surv() object
    names(data) <- varnames

    ## Any remaining matrices not of class 'rms' must be given class rms
    ## so that as.data.frame will not split up their variables
    ism <- sapply(data, ismat)
    if(any(ism))
      for(i in which(ism))
        class(data[[i]]) <- unique(c('rms', class(data[[i]])))

    ## Since subsetting was completed earlier, now drop unused factor levels
    ## NOTE: strat() variables are also factors; don't drop their attributes
    isf <- sapply(data, is.factor)
    if(any(isf))
      for(i in which(isf)) {
        di <- data[[i]]
        at <- attributes(di)
        di <- di[, drop=TRUE]
        if(length(at$assume.code) && at$assume.code == 8) {
          at$levels   <- at$parms <- levels(di)
          at$colnames <- paste0(at$name, '=', levels(di)[-1])
          attributes(di) <- at[c('class', 'name', 'label', 'assume',
                                 'assume.code', 'parms', 'nonlinear', 'tex',
                                 'colnames','levels')]
          data[[i]] <- di
        }
      }

    ## If any variables are less than the maximum length, these must
    ## have come from the parent environment and did not have subset applied
    len <- sapply(data, NROW)
    if(min(len) != max(len)) {
      if(! length(subset))
        stop('program logic error: variables vary in length but subset= was not given')
      for(i in which(len > min(len))) {
        x <- data[[i]]
        data[[i]] <- if(is.matrix(x)) x[subset,,drop=FALSE] else x[subset]
      }
      len <- sapply(data, NROW)
      if(min(len) != max(len)) stop('program logic error in variable lengths')
    }
    # row.names added 2023-01-25
    data        <- as.data.frame(data, check.names=FALSE, row.names=rnames)
    data        <- na.action(data)
    nac         <- attr(data, 'na.action')
    attr(data, 'na.action') <- nac
    data
  }
  dat  <- processdata(formula, data)
  if(length(formula2)) {
    omit <- attr(dat, 'na.action')$omit
    if(length(omit)) data <- data[-omit, , drop=FALSE]
    dat2 <- processdata(formula2, data)
    attr(dat, 'data2') <- dat2
    }
  dat
  }

## Handle spline and other variables with rms class
as.data.frame.rms <- function(x, row.names = NULL, optional = FALSE, ...) {
  nrows <- NROW(x)
  row.names <- if(optional) character(nrows) else as.character(1:nrows)
  value <- list(x)
  if(! optional) names(value) <- deparse(substitute(x))[[1]]
  structure(value, row.names=row.names, class='data.frame')
}
