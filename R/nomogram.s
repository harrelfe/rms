nomogram <-
  function(fit, ..., adj.to, 
           lp=TRUE, lp.at=NULL,
           fun=NULL, fun.at=NULL, fun.lp.at=NULL, funlabel="Predicted Value",
           interact=NULL, kint=NULL,
           conf.int=FALSE, 
           conf.lp=c("representative", "all", "none"),
           est.all=TRUE, abbrev=FALSE, minlength=4, maxscale=100, nint=10, 
           vnames=c("labels","names"),
           varname.label=TRUE, varname.label.sep="=",
           omit=NULL, verbose=FALSE)
{	

  conf.lp <- match.arg(conf.lp)
  vnames  <- match.arg(vnames)

  Format <- function(x)
    { # like format but does individually
      f <- character(l <- length(x))
      for(i in 1:l) f[i] <- format(x[i])
      f
    }

  abb <- (is.logical(abbrev) && abbrev) || is.character(abbrev)
  if(is.logical(conf.int) && conf.int) conf.int <- c(.7,.9)

  se <- any(conf.int > 0)

  nfun <- if(!length(fun)) 0 else if(is.list(fun)) length(fun) else 1
  if(nfun>1 && length(funlabel) == 1) funlabel <- rep(funlabel, nfun)
  if(nfun>0 && is.list(fun) && length(names(fun))) funlabel <- names(fun)

  if(length(fun.at) && !is.list(fun.at)) 
    fun.at <- rep(list(fun.at), nfun)
  if(length(fun.lp.at) && !is.list(fun.lp.at))
    fun.lp.at <- rep(list(fun.lp.at), nfun)

  at <- fit$Design
  assume <- at$assume.code
  if(any(assume == 10))
    warning("does not currently work with matrix factors in model")
  name  <- at$name
  names(assume) <- name
  parms <- at$parms
  label <- if(vnames == "labels") at$label else name
  if(any(d <- duplicated(name)))
    stop(paste("duplicated variable names:",
               paste(name[d],collapse=" ")))
  label <- name
  if(vnames == "labels") {
    label <- at$label
    if(any(d <- duplicated(label)))
      stop(paste("duplicated variable labels:",
                 paste(label[d],collapse=" ")))
  }
  
  ia    <- at$interactions
  
  factors <- rmsArgs(substitute(list(...)))
  nf <- length(factors)

  which <- if(est.all) (1:length(assume))[assume != 8]
  else 
    (1:length(assume))[assume != 8 & assume != 9]
  if(nf > 0) {
    jw <- charmatch(names(factors), name, 0)
    if(any(jw == 0)) stop(paste("factor name(s) not in the design:",
             paste(names(factors)[jw == 0], collapse=" ")))
    if(!est.all) which <- jw
  }

  Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
  values <- Limval$values
  lims   <- Limval$limits[c(6, 2, 7),, drop=FALSE]

  ## Keep character variables intact
  lims <- unclass(lims)
  for(i in 1:length(lims))
    if(is.factor(lims[[i]])) lims[[i]] <- as.character(lims[[i]])
  attr(lims, 'class') <- 'data.frame'  # so can subscript later

  ## Find underlying categorical variables
  ucat <- rep(FALSE, length(assume))
  names(ucat) <- name
  for(i in (1:length(assume))[assume != 5 & assume < 8]) {
    ucat[i] <- !is.null(V <- values[[name[i]]]) # did add && is.character(V)
    if(ucat[i]) parms[[name[i]]] <- V
  }

  discrete <- assume == 5 | assume == 8 | ucat
  names(discrete) <- name

  ## Number of non-slopes:
  nrp <- num.intercepts(fit, 'coef')
  ir <- fit$interceptRef
  if(!length(ir))   ir   <- 1
  if(!length(kint)) kint <- ir
  
  Intercept <- if(nrp > 0) fit$coefficients[kint]
  else 
    if(length(fit$center)) (- fit$center ) else 0
  intercept.offset <- fit$coefficients[kint] - fit$coefficients[ir]

  settings <- list()
  for(i in which[assume[which] < 9]) {
    ni <- name[i]
    z <- factors[[ni]]
    lz <- length(z)
    if(lz < 2) settings[[ni]] <- value.chk(at, i, NA, -nint, Limval, 
                                           type.range="full") else
    if(lz > 0 && any(is.na(z)))
      stop("may not specify NA as a variable value")
    if(lz == 1) lims[2,i] <- z else if(lz > 1) {
      settings[[ni]] <- z
      if(is.null(lims[[ni]]) || is.na(lims[2, ni])) {
        lims[[ni]] <- c(NA, z[1], NA)
        warning(paste("adjustment values for ",ni,
                      " not defined in datadist; taken to be first value specified (", 
                      z[1], ")" ,sep=""))
      }
    }
  }

  adj <- lims[2,, drop=FALSE]
  if(!missing(adj.to))
    for(nn in names(adj.to)) adj[[nn]] <- adj.to[[nn]]
  isna <- sapply(adj, is.na)
  if(any(isna))
    stop(
      paste("adjustment values not defined here or with datadist for",
            paste(name[assume != 9][isna],collapse=" ")))
  
  num.lines <- 0

  entities <- 0
  main.space.used <- ia.space.used <- 0

  set  <- list()
  nset <- character(0)
  iset <- 0

  start <- len <- NULL
  end <- 0

  ## Sort to do continuous factors first if any interactions present
  
  main.effects <- which[assume[which] < 8]
  ## this logic not handle strata w/intera.
  if(any(assume == 9))
    main.effects <-
      main.effects[order(10 * discrete[main.effects] +
                         (name[main.effects] %in% names(interact)))]
  
  ## For each predictor, get vector of predictor numbers directly or
  ## indirectly associated with it
  rel <- related.predictors(at)   # Function in rmsMisc.s

  already.done <- structure(rep(FALSE,length(name)), names=name)
  for(i in main.effects) {
    nam <- name[i]
    if(already.done[nam] || (nam %in% omit)) next
    r <- if(length(rel[[nam]])) sort(rel[[nam]]) else NULL
    if(length(r) == 0) { #main effect not contained in any interactions
      num.lines  <- num.lines + 1
      main.space.used <- main.space.used + 1
      entities   <- entities + 1
      x <- list()
      x[[nam]] <- settings[[nam]]
      iset <- iset + 1
      attr(x,'info') <- list(nfun=nfun, predictor=nam, effect.name=nam,
                             type='main')
      set[[iset]] <- x
      nset <- c(nset, label[i])
      
      start <- c(start, end + 1)
      n <- length(settings[[nam]])
      len <- c(len, n)
      end <- end + n
    } else {
      namo <- name[r]
      s <- !(name[r] %in% names(interact))
      if(any(s)) {
        if(!length(interact)) interact <- list()
        for(j in r[s]) {
          nj <- name[j]
          if(discrete[j]) interact[[nj]] <-  parms[[nj]]
        }
        s <- !(name[r] %in% names(interact))
      }
      if(any(s)) stop(paste("factors not defined in interact=list(...):",
                            paste(name[r[s]], collapse=",")))
      combo <- expand.grid(interact[namo]) #list[vector] gets sublist
      class(combo) <- NULL
      ## so combo[[n]] <- as.character will really work
      acombo <- combo
      if(abb)
        for(n in if(is.character(abbrev)) abbrev else names(acombo)) {
          if(discrete[n]) { 
            acombo[[n]] <-
              abbreviate(parms[[n]],
                         minlength=if(minlength == 1) 4
                                   else minlength)[combo[[n]]]
            ## lucky that abbreviate function names its result
          }
        }
      for(n in names(combo)) if(is.factor(combo[[n]])) {
        combo[[n]] <- as.character(combo[[n]])
        ## so row insertion will work xadj
        acombo[[n]] <- as.character(acombo[[n]]) #so format() will work
      }
      entities <- entities + 1
      already.done[namo] <- TRUE
      for(k in 1:length(combo[[1]])) {
        num.lines <- num.lines + 1
        if(k == 1) main.space.used <- main.space.used + 1
        else ia.space.used <- ia.space.used + 1
        x <- list()
        x[[nam]] <- settings[[nam]]   #store fastest first
        for(nm in namo) x[[nm]] <- combo[[nm]][k]
        iset <- iset + 1
        set.name <- paste(nam, " (", sep="")
        for(j in 1:length(acombo)) {
          set.name <-
            paste(set.name, 
                  if(varname.label) paste(namo[j], varname.label.sep,
                                          sep="") else "",
                  format(acombo[[j]][k]), sep="")
          if(j < length(acombo)) set.name <- paste(set.name," ",sep="")
        }
        set.name <- paste(set.name, ")", sep="")
        ## Make list of all terms needing inclusion in calculation
        ## Include interation term names  - interactions.containing in rmsMisc.s
        ia.names <- NULL
        for(j in r) ia.names <-
          c(ia.names, name[interactions.containing(at, j)])
        ia.names <- unique(ia.names)
        attr(x,'info') <-
          list(predictor=nam,
               effect.name=c(nam, namo[assume[namo] != 8], ia.names), 
               type=if(k == 1) "first" else "continuation")
        set[[iset]] <- x
        nset <- c(nset, set.name)
        ## Don't include strata main effects
        start <- c(start, end + 1)
        n <- length(settings[[nam]])
        len <- c(len, n)
        end <- end + n
      }
    }    
  }
  xadj <- unclass(rms.levels(adj, at))
  for(k in 1:length(xadj)) xadj[[k]] <- rep(xadj[[k]], sum(len))
  
  j <- 0
  for(S in set) {
    j <- j + 1
    ns  <- names(S)
    nam <- names(S)
    for(k in 1:length(nam))
      xadj[[nam[k]]][start[j] : (start[j] + len[j]-1)] <- S[[k]]
  }
  xadj <- structure(xadj, class='data.frame',
                    row.names=as.character(1 : sum(len)))

  xx <- predictrms(fit, newdata=xadj, type="terms",
                   center.terms=FALSE, se.fit=FALSE, kint=kint)
  if(any(is.infinite(xx)))
    stop("variable limits and transformations are such that an infinite axis value has resulted.\nRe-run specifying your own limits to variables.")
  
  if(se) xse <- predictrms(fit, newdata=xadj, se.fit=TRUE, 
                           kint=kint)

  R <- matrix(NA, nrow=2, ncol=length(main.effects),
              dimnames=list(NULL,name[main.effects]))
  R[1,] <-  1e30
  R[2,] <- -1e30
  ## R <- apply(xx, 2, range)  - does not work since some effects are for
  ## variable combinations that were never used in constructing axes

  for(i in 1:num.lines) {
    is <- start[i]; ie <- is + len[i]-1
    s <- set[[i]]
    setinfo <- attr(s, 'info')
    nam <- setinfo$effect.name
    xt <- xx[is : ie, nam]
    if(length(nam) > 1) xt <- apply(xt, 1, sum)  # add all terms involved
    set[[i]]$Xbeta <- xt
    r <- range(xt)
    pname <- setinfo$predictor
    R[1,pname] <- min(R[1,pname], r[1])
    R[2,pname] <- max(R[2,pname], r[2])
    if(se) {
      set[[i]]$Xbeta.whole <-
        xse$linear.predictors[is:ie] #note-has right interc.
      set[[i]]$se.fit      <- xse$se.fit[is:ie]
    }
  }

  R <- R[,R[1,] < 1e30, drop=FALSE]
  sc <- maxscale / max(R[2,] - R[1,])
  Intercept <- Intercept + sum(R[1,])
  
###if(missing(naxes)) naxes <- 
###  if(total.sep.page) max(space.used + 1, nfun + lp + 1) else
###                     space.used + 1 + nfun + lp + 1
  
  i <- 0
  names(set) <- nset
  ns <- names(set)
  Abbrev <- list()
  qualForce <- character()
  for(S in set) {
    i <- i + 1
    setinfo <- attr(S,'info')
    type <- setinfo$type
    x <- S[[1]]
    nam <- names(S)[1]  #stored with fastest first
    fx <- if(is.character(x)) x else 
    sedit(Format(x)," ","") #axis not like bl   - was translate()
    if(abb && discrete[nam] && (is.logical(abbrev) || nam %in% abbrev)) {
      old.text <- fx
      fx <- if(abb && minlength == 1)letters[1:length(fx)] else 
      abbreviate(fx, minlength=minlength)
      Abbrev[[nam]] <- list(abbrev=fx, full=old.text)
    }
    
    j <- match(nam, name, 0)
    if(any(j == 0)) stop("program logic error 1")
    is <- start[i]
    ie <- is + len[i] - 1
    xt <- (S$Xbeta - R[1,nam]) * sc
    set[[i]]$points <- xt
    ## Find flat pieces and combine their labels
    r <- rle(xt)
    if(any(r$length > 1)) {
      is <- 1
      for(j in r$length) {
        ie <- is + j - 1
        if(j > 1) {
          fx[ie] <- if(discrete[nam] || ie < length(xt))
            paste(fx[is], "-", fx[ie], sep="") else
          paste(fx[is], '+', sep='')
          
          fx[is:(ie-1)] <- ""
          xt[is:(ie-1)] <- NA
        }
        is <- ie + 1
      }
      fx <- fx[!is.na(xt)]
      xt <- xt[!is.na(xt)]
    }
  }
  if(!length(lp.at)) {
    xb <- fit$linear.predictors
    if(!length(xb)) xb <- fit$fitted.values
    if(!length(xb)) xb <- fit$fitted
    if(!length(xb))
      stop("lp.at not given and fit did not store linear.predictors or fitted.values")
    if(nrp > 1) xb <- xb + intercept.offset
    lp.at <- pretty(range(xb), n=nint)
  }
  
  sum.max <- if(entities == 1) maxscale
  else max(maxscale, sc * max(lp.at - Intercept))
  x <- pretty(c(0, sum.max), n=nint)
  
  new.max <- max(x)
  
  iset <- iset + 1
  nset <- c(nset, 'total.points')
  set[[iset]] <- list(x=x)
  
  if(lp) {
    x2 <- seq(lp.at[1], max(lp.at), by=(lp.at[2] - lp.at[1]) / 2)
    scaled.x  <- (lp.at - Intercept) * sc
    iset <- iset + 1
    nset <- c(nset, 'lp')
    if(se && conf.lp != 'none') {
      xxb <- NULL
      xse <- NULL
      for(S in set) {
        xxb <- c(xxb, S$Xbeta.whole)
        xse <- c(xse, S$se.fit)
      }
      i <- order(xxb)
      if(length(xxb)<16 | conf.lp == "representative") 
        {nlev <- 4; w <- 1} else {nlev <- 8; w <- 2}
      if(conf.lp == "representative") {
        deciles   <- cut2(xxb[i], g=10)
        mean.xxb  <- tapply(xxb[i], deciles, mean)
        median.se <- tapply(xse[i], deciles, median)
        xc  <- (mean.xxb - Intercept)*sc
        sec <- sc*median.se
      }
      else {
        xc  <- (xxb[i]-Intercept)*sc
        sec <- sc*xse[i]
      }
      set[[iset]] <- list(x=scaled.x, x.real=lp.at,
                          conf=list(x=xc, se=sec, w=w, nlev=nlev))
    }
    else set[[iset]] <- list(x=scaled.x, x.real=lp.at)
  }
  
  if(nfun > 0) {
    if(!is.list(fun)) fun <- list(fun)
    i <- 0
    for(func in fun) {
      i <- i + 1
      ## Now get good approximation to inverse of fun evaluated at fat
      ## Unless inverse function given explicitly
      if(!missing(fun.lp.at)) {
        xseq <- fun.lp.at[[i]]
        fat  <- func(xseq)
        w    <- xseq
      } else {
        if(missing(fun.at)) fat <- pretty(func(range(lp.at)), n=nint)
        else fat <- fun.at[[i]]
        if(verbose) {
          cat('Function',i,'values at which to place tick marks:\n')
          print(fat)
        }
        xseq <- seq(min(lp.at), max(lp.at), length=1000)
        fu   <- func(xseq)
        s    <- !is.na(fu)
        w    <- approx(fu[s], xseq[s], fat, ties=mean)$y
        if(verbose) {
          cat('Estimated inverse function values (lp):\n')
          print(w)
        }
      }
      s <- !(is.na(w) | is.na(fat))
      w <- w[s]
      fat <- fat[s]
      fat.orig <- fat
      fat <- if(is.factor(fat)) as.character(fat) else Format(fat)
      scaled <- (w - Intercept) * sc
      iset <- iset + 1
      nset <- c(nset, funlabel[i])
      set[[iset]] <- list(x=scaled, x.real=fat.orig, fat=fat, which=s)
    }
  }
  names(set) <- nset
  attr(set, 'info') <- list(fun=fun, lp=lp, lp.at=lp.at,
                            discrete=discrete, funlabel=funlabel,
                            fun.at=fun.at, fun.lp.at=fun.lp.at,
                            Abbrev=Abbrev,  minlength=minlength,
                            conf.int=conf.int,
                            R=R, sc=sc, maxscale=maxscale,
                            Intercept=Intercept, nint=nint,
                            space.used=c(main=main.space.used, ia=ia.space.used))
  class(set) <- "nomogram"
  set
}

print.nomogram <- function(x, dec=0, ...)
{
  obj <- x
  w <- diff(range(obj$lp$x)) / diff(range(obj$lp$x.real))
  cat('Points per unit of linear predictor:', format(w),
	  '\nLinear predictor units per point   :', format(1 / w), '\n\n')
  
  fun <- FALSE
  for(x in names(obj)) {
    k <- x == 'total.points' || x == 'lp' || x == 'abbrev'
    if(k) { fun <- TRUE; next }
    y <- obj[[x]]
    if(fun) {
      z <- cbind(round(y[[1]],dec), y$x.real)
      dimnames(z) <- list(rep('',nrow(z)), c('Total Points',x))
    } else {
      z <- cbind(format(y[[1]]), format(round(y$points,dec)))
      dimnames(z)  <- list(rep('',length(y$points)), c(x, 'Points'))
      ## didn't use data.frame since wanted blank row names
    }
    cat('\n')
    print(z, quote=FALSE)
    cat('\n')
  }
  invisible()
}
