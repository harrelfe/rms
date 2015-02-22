plot.Predict <-
  function(x, formula, groups=NULL, cond=NULL, varypred=FALSE, subset,
           xlim, ylim, xlab, ylab,
           data=NULL, subdata, anova=NULL, pval=FALSE, cex.anova=.85,
           col.fill=gray(seq(.825, .55, length=5)),
           adj.subtitle, cex.adj, cex.axis, perim=NULL,
           digits=4, nlevels=3, nlines=FALSE, addpanel,
           scat1d.opts=list(frac=0.025, lwd=0.3), type=NULL,
           yscale=NULL, scaletrans=function(z) z, ...)
{
  if(varypred) {
    x$.predictor. <- x$.set.
    x$.set. <- NULL
  }
  predpres <- length(x$.predictor.)
  if(missing(addpanel)) addpanel <- function(...) {}

  subdatapres <- !missing(subdata)
  if(subdatapres) subdata <- substitute(subdata)
  
  doscat1d <- function(x, y, col) {
    so <- scat1d.opts
    if(!length(so$col)) so$col <- col
    do.call('scat1d', c(list(x=x, y=y), so, grid=TRUE))
  }
  
  info     <- attr(x, 'info')
  at       <- info$Design
  label    <- at$label
  units    <- at$units
  values   <- info$values
  adjust   <- info$adjust
  yunits   <- info$yunits
  varying  <- info$varying
  conf.int <- info$conf.int
  dotlist  <- list(...)

  gname <- groups
  if(length(gname)) {
    if(length(gname) > 1 || !is.character(gname) ||
       gname %nin% names(x))
      stop('groups must be a single predictor name')
  }
  if(length(cond)) {
    if(length(cond) > 1 || !is.character(cond) ||
       cond %nin% names(x))
      stop('cond must be a single predictor name')
  }
  if(missing(ylab))    ylab     <- info$ylabPlotmath
  if(!length(x$lower)) conf.int <- FALSE
  
  if(missing(ylim))
    ylim <- range(pretty(
             if(conf.int) c(x$yhat, x$lower, x$upper)
              else x$yhat), na.rm=TRUE)

  if(missing(adj.subtitle)) adj.subtitle <- length(adjust) > 0
  sub <- if(adj.subtitle && length(adjust)==1)
    paste('Adjusted to:', adjust, sep='') else NULL
  cex <- par('cex')
  if(missing(cex.adj)) cex.adj <- .75*cex
  if(length(sub)) sub <- list(sub, cex=cex.adj, just=c('center','bottom'))
  
  subset <- if(missing(subset)) TRUE
  else
    eval(substitute(subset),x)

  oldopt <- options(digits=digits)
  on.exit(options(oldopt))

  if(length(anova)) {
    stat   <- plotmathAnova(anova, pval)
    tanova <- function(name, x, y)
      annotateAnova(name, stat, x, y, cex=cex.anova)
  }
  else tanova <- function(...) {}
    
  if(predpres) {
    if(! missing(formula))
      stop('formula may not be given when multiple predictors were varied separately')
    p  <- as.factor(x$.predictor.)
    xp <- rep(NA, length(p))
    levs <- at <- labels <- limits <- list()
    lp <- levels(p)
    np <- length(lp)
    groups <- gr <- if(length(gname)) as.factor(x[[gname]])
    cond   <- co <- if(length(cond))  as.factor(x[[cond]])
    
    perhapsAbbrev <- function(k) {
      len <- nchar(k)
      if(sum(len) > 30 | max(len) > 12)
        abbreviate(k, minlength=max(3, round(17 / length(k))))
      else
        k
    }
    for(w in lp) {
      i  <- p == w
      z  <- x[[w]]
      l  <- levels(z)
      ll <- length(l)
      levs[[w]]   <- if(ll) l else character(0)
      xp[i]       <- as.numeric(z[i])
      if(length(groups)) gr[i] <- groups[i]
      if(length(cond))   co[i] <- cond[i]
      at[[w]]     <- if(ll) 1 : ll else pretty(z[i])
      labels[[w]] <- if(ll) perhapsAbbrev(l) else format(at[[w]])
      limits[[w]] <- if(ll) c(2 / 3, ll + 1 / 3) else range(z[i])
    }
    if(length(cond)) {
      nuc <- length(levels(cond))
      at <- at[rep(seq(1, length(at)), each=nuc)]
      labels <- labels[rep(seq(1, length(labels)), each=nuc)]
      limits <- limits[rep(seq(1, length(limits)), each=nuc)]
      levs <- levs[rep(seq(1, length(levs)), each=nuc)]
      
      formula <- if(!conf.int) x$yhat ~ xp | cond*p else
      Cbind(x$yhat, x$lower, x$upper) ~ xp | cond*p
    }
    else {
      formula <- if(!conf.int) x$yhat ~ xp | p else
      Cbind(x$yhat, x$lower, x$upper) ~ xp | p
    }
    
    pan <- function(x, y, ...) {
      pn <- lattice::panel.number()
      lev    <- levs[[pn]]
      col <- lattice::trellis.par.get('superpose.line')$col
      if(!length(lev) && length(unique(x[!is.na(x)])) > nlevels)
        { # continuous x
          yy <- y
          if(length(perim)) {
            j <- perim(x, NULL)
            yy[j] <- NA
            if(length(attr(yy, 'other'))) attr(yy, 'other')[j, ] <- NA
          }
          panel.xYplot(x, yy, ...)
          tanova(names(levs)[pn], x, yy)
          if(length(data) && length(xd <- data[[names(levs)[pn]]])) {
            xd <- xd[!is.na(xd)]
            doscat1d(xd, approx(x, y, xout=xd, rule=2, ties=mean)$y,
                     col=col[1])
          }
        }
      else
        { # discrete x
          lattice::panel.points(x, y, pch=19)
          yoth <- attr(y, 'other')
          yo   <- length(yoth)
          if(yo) for(u in unique(x))
            lattice::llines(c(u, u), yoth[x==u, ])
          tanova(names(levs)[pn],
                 if(yo) c(x,         x,         x) else x,
                 if(yo) c(y, yoth[, 1], yoth[, 2]) else y)
        }
      addpanel(x, y, ...)
    }
    scales <- list(x=list(relation='free', limits=limits,
                     at=at, labels=labels))
    if(!missing(cex.axis)) scales$x$cex <- cex.axis
    if(length(yscale)) scales$y <- yscale
    r <- list(formula=formula, groups=gr, subset=subset,
              type=if(length(type))type else 'l',
              method=if(conf.int & (!length(type) || type != 'p'))
              'filled bands' else 'bars',
              col.fill=col.fill,
              xlab='', ylab=ylab, ylim=ylim,
              panel=pan, scales=scaletrans(scales),
              between=list(x=.5))
    if(length(dotlist)) r <- c(r, dotlist)
    if(length(sub   )) r$sub    <- sub
  }
  else
    { # predictor not listed
      v <- character(0)
      bar <- ''
      f <- if(!missing(formula)) gsub(' .*','',as.character(formula)[2])
      else varying[1]
      iv <- var.inner(as.formula(paste('~', f)))
      if(missing(xlab)) xlab <- labelPlotmath(label[iv], units[iv])
      if(missing(formula)) {
        xvar <- varying[1]
        ## change formula like ~x|foo to x
        if(length(varying) > 1) {
          v <- varying[-1]
          if(length(gname)) {
            groups <- x[[gname]]
            v <- setdiff(v, gname)
          }
          else {
            nu <- sapply(x[v], function(u) length(unique(u)))
            if(!predpres && any(nu <= nlevels)) {
              i <- which.min(nu)
              gname <- v[i]
              groups <- x[[gname]]
              v <- setdiff(v, gname)
            }
          }
          if(length(v)) {
            bar <- paste(v, collapse='*')
            for(z in v) {
              if(z %in% names(x) && is.numeric(x[[z]])) {
                x[[z]] <- factor(format(x[[z]]))
                levels(x[[z]]) <- paste(z,':',levels(x[[z]]), sep='')
              }
            }
          }
        } 
        xv <- x[[xvar]]
        xdiscrete <- is.factor(xv) || is.character(xv) ||
          length(unique(xv[!is.na(xv)])) <= nlevels
        if(xdiscrete) {
          f <- paste(xvar, if(conf.int) 'Cbind(yhat,lower,upper)'
          else 'yhat', sep='~')
          if(bar != '') f <- paste(f, bar, sep='|')
          formula <- eval(parse(text=f))
          if(length(v)) for(z in v) {
            if(z %in% names(x) && is.numeric(x[[z]])) {
              x[[z]] <- factor(format(x[[z]]))
              levels(x[[z]]) <- paste(z,':',levels(x[[z]]), sep='')
            }
          }
          r <- Dotplot(formula, groups=groups, subset=subset,
                       xlim=ylim, xlab=ylab, ylab=xlab,
                       sub=sub, data=x, between=list(x=.5), ...)
          return(r)
        }
        if(bar != '') f <- paste(f, '|', bar)
      }
      else
        { # formula given
          f <- as.character(formula)[2]
          xvar <- gsub(' .*', '', f)
          if(length(grep('\\|', f))) {
            g <- gsub(' ', '', f)
            g <- gsub('.*\\|', '', g)
            v <- strsplit(g, '\\*')[[1]]
          }
          if(length(v)) for(z in v) {
            if(z %in% names(x) && is.numeric(x[[z]])) {
              x[[z]] <- factor(format(x[[z]]))
              levels(x[[z]]) <- paste(z,':',levels(x[[z]]), sep='')
            }
          }
        }
      f <- paste(if(conf.int) 'Cbind(yhat,lower,upper)' else 'yhat',
                 f, sep='~')
      formula <- eval(parse(text=f))
      
      xv <- x[[xvar]]
      xscale <- NULL
      xdiscrete <- (is.factor(xv) || is.character(xv)) && nlines
      if(xdiscrete) {
        xv <- as.factor(xv)
        xlev <- levels(xv)
        xscale <- list(x=list(at=1:length(xlev), labels=xlev))
        if(!missing(cex.axis)) xscale$x$cex <- cex.axis
        x[[xvar]] <- as.integer(xv)
      }

      ## Continuing: no predpres case
      pan <- function(x, y, groups=NULL, subscripts, ...) {
        ogroups <- groups
        if(length(groups)) groups <- groups[subscripts]
        yy <- y
        if(length(perim)) {
          if(! length(groups)) {
            j <- ! perim(x, NULL)
            yy[j] <- NA
            if(length(attr(yy, 'other'))) attr(yy, 'other')[j, ] <- NA
          }
          else { ## perim and groups specified
            for(w in if(is.factor(groups)) levels(groups) else unique(groups)) {
              i  <- which(groups == w)
              j <- ! perim(x[i], w)
              yy[i[j]] <- NA
              if(length(attr(yy, 'other')))
                attr(yy, 'other')[i[j], ] <- NA
            }
          }
        }
        panel.xYplot(x, yy, groups=ogroups, subscripts=subscripts, ...)
        tanova(xvar, x, yy)
        col <- lattice::trellis.par.get('superpose.line')$col
        
        xd <- data[[xvar]]
        use <- TRUE
        if(length(xd) && subdatapres) {
          use <- eval(subdata, data)
          if(length(use) != nrow(data))
            stop('subdata must evaluate to a length of nrow(data)')
        }
        
        if(length(groups) && length(gd <- data[[gname]]) &&
           length(xd)) {
          g <- groups[subscripts]
          j <- 0
          for(w in levels(g)) {
            j <- j + 1
            z  <- g==w
            xg <- x[z]
            yg <- y[z]
            x1 <- xd[use & gd==w]
            x1 <- x1[!is.na(x1)]
            doscat1d(x1, approx(xg, yg, xout=x1, rule=2, ties=mean)$y,
                     col=col[j])
          }
        }
        else if(length(xd)) {
          xd <- xd[use & !is.na(xd)]
          doscat1d(xd, approx(x, y, xout=xd, rule=2, ties=mean)$y,
                   col=col[1])
        }
        addpanel(x, y, groups=NULL, subscripts=subscripts, ...)
      }
      
      r <- list(formula=formula, data=x, subset=subset,
                type=if(length(type)) type else if(xdiscrete) 'b' else 'l',
                method=if(conf.int & (!length(type) || type!='p'))
                       'filled bands' else 'bars',
                col.fill=col.fill,
                xlab=xlab, ylab=ylab, ylim=ylim, panel=pan,
                between=list(x=.5))
      scales <- NULL
      if(length(xscale)) scales <- xscale
      if(length(yscale)) scales$y <- yscale
      r$scales <- scaletrans(scales)
      if(!missing(xlim)) r$xlim   <- xlim
      if(!conf.int)      r$method <- NULL
      if(length(gname))  r$groups <- x[[gname]]
      if(length(sub))    r$sub    <- sub
      if(length(dotlist)) r <- c(r, dotlist)
    }
  do.call('xYplot', r)
}


pantext <- function(object, x, y, cex=.5, adj=0,
                    fontfamily='Courier', lattice=TRUE)
{
  k <- paste(capture.output(object), collapse='\n')
  fam <- fontfamily
  if(lattice) {
    z <- 
      function(x, y, ..., xx, yy, text, cex, adj, family)
        lattice::ltext(xx, yy, text, cex=cex, adj=adj, fontfamily=family)
    formals(z) <-
      eval(substitute(alist(x=, y=, ...=, xx=xx, yy=yy, text=k,
                            cex=cex, adj=adj, family=fam),
                      list(xx=x, yy=y, k=k, cex=cex, adj=adj,
                           fam=fam)))
  }
  else {
    z <-
      function(x, y, text, cex, adj, family, ...)
        text(x, y, text, adj=adj, cex=cex, family=family, ...)
    formals(z) <-
      eval(substitute(alist(x=x, y=y, text=k,
                            cex=cex, adj=adj, family=fam, ...=),
                      list(x=x, y=y, k=k, cex=cex, adj=adj,
                           fam=fam)))
  }
  z
}

plotmathAnova <- function(anova, pval) {
  vi <- attr(anova, 'vinfo')
  aname  <- sapply(vi, function(x) paste(x$name, collapse=','))
  atype  <- sapply(vi, function(x) x$type)
  wanova <- atype %in% c('main effect', 'combined effect')
  test   <- if('F' %in% colnames(anova)) 'F' else 'Chi-Square'
  stat   <- round(anova[wanova, test], 1)
  pstat  <- anova[wanova, 'P']
  dof    <- anova[wanova, 'd.f.']
  stat <- if(test == 'Chi-Square')
    paste('chi[', dof, ']^2 == ',  stat, sep='')
  else
    paste('F[paste(', dof, ',",",', anova['ERROR', 'd.f.'], ')] == ',
          stat, sep='')
  if(pval) {
    pval <- formatNP(pstat, digits=3, pvalue=TRUE)
    pval <- ifelse(grepl('<', pval), paste('P', pval, sep=''),
                   paste('P==', pval, sep=''))
    stat <- paste(stat, pval, sep='~~')
  }
  names(stat) <- aname[wanova]
  stat
}

## stat is result of plotmathAnova
## xlim and ylim must be specified if ggplot=TRUE
annotateAnova <- function(name, stat, x, y, ggplot=FALSE,
                          xlim=NULL, ylim=NULL, cex, size=4, flip=FALSE,
                          empty=FALSE, dataOnly=FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if(flip) {
    yorig <- y
    y     <- x
    x     <- yorig
    ylimorig <- ylim
    ylim     <- xlim
    xlim     <- ylimorig
  }
    
  ## size is for ggplot2 only; is in mm
  ## See if an area is available near the top or bottom of the
  ## current panel
  if(! ggplot) {
    cpl <- lattice::current.panel.limits()
    xlim <- cpl$xlim
    ylim <- cpl$ylim
  }
  else if(! length(xlim) || ! length(ylim))
    stop('xlim and ylim must be given if ggplot=TRUE')
  dy   <- diff(ylim)
  if(!empty && !any(y > ylim[2] - dy / 7)) {
    z <- list(x = mean(xlim), y = ylim[2] - .025 * dy)
    adj <- c(.5, 1)
  }
  else if(! empty && !any(y < ylim[1] + dy / 7)) {
    z <- list(x = mean(xlim), y = ylim[1] + .025 * dy)
    adj <- c(.5, 0)
  }
  else {
    z <- if(! length(xlim) || ! length(ylim)) 
      largest.empty(x, y, grid=TRUE, method='exhaustive')
    else
      largest.empty(x, y, grid=TRUE, method='exhaustive', xlim=xlim, ylim=ylim)
    adj <- c(if(z$x > mean(xlim)) 1 else .5,
             if(z$y > mean(ylim)) 1 else 0)
  }
  if(flip) {
    zyorig <- z$y
    z$y <- z$x
    z$x <- zyorig
    adj <- rev(adj)
  }
  ## parse=TRUE: treat stat[name] as an expression
  if(dataOnly) return(list(x=z$x, y=z$y, label=stat[name],
                           hjust=adj[1], vjust=adj[2]))
  if(ggplot) annotate('text', x=z$x, y=z$y, label=stat[name], parse=TRUE,
                       size=size, hjust=adj[1], vjust=adj[2])
  else lattice::ltext(z$x, z$y, parse(text=stat[name]), cex=cex, adj=adj)
}
