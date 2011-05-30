plot.Predict <-
  function(x, formula, groups=NULL, cond=NULL, varypred=FALSE, subset,
           xlim, ylim, xlab, ylab,
           data=NULL, col.fill=gray(seq(.95, .75, length=5)),
           adj.subtitle, cex.adj, perim=NULL,
           digits=4, nlevels=3, nlines=FALSE, addpanel,
           scat1d.opts=list(frac=0.025, lwd=0.3), ...)
{
  require(lattice)
  if(varypred)
    {
      x$.predictor. <- x$.set.
      x$.set. <- NULL
    }
  predpres <- length(x$.predictor.)
  if(missing(addpanel)) addpanel <- function(...) {}

  doscat1d <- function(x, y, col)
    {
      so <- scat1d.opts
      if(!length(so$col)) so$col <- col
      do.call('scat1d', c(list(x=x, y=y), so, grid=TRUE))
    }

  info   <- attr(x, 'info')
  at     <- info$Design
  label  <- at$label
  units  <- at$units
  values <- info$values
  adjust <- info$adjust
  yunits <- info$yunits
  varying<- info$varying
  conf.int <- info$conf.int
  dotlist <- list(...)

  gname <- groups
  if(length(gname))
    {
      if(length(gname) > 1 || !is.character(gname) ||
         gname %nin% names(x))
        stop('groups must be a single predictor name')
    }
  if(length(cond))
    {
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

  if(predpres)
    {
      if(!missing(formula))
        stop('formula may not be given when multiple predictors were varied separately')
      p  <- as.factor(x$.predictor.)
      xp <- rep(NA, length(p))
      levs <- at <- labels <- limits <- list()
      lp <- levels(p)
      np <- length(lp)
      groups <- gr <- if(length(gname)) as.factor(x[[gname]])
      cond <- co <- if(length(cond)) as.factor(x[[cond]])

      perhapsAbbrev <- function(k)
        {
          len <- nchar(k)
          if(sum(len) > 30 | max(len) > 12)
            abbreviate(k, minlength=max(3,round(17/length(k))))
          else
            k
        }
      for(w in lp)
        {
          i  <- p==w
          z  <- x[[w]]
          l  <- levels(z)
          ll <- length(l)
          levs[[w]]   <- if(ll) l else character(0)
          xp[i]       <- as.numeric(z[i])
          if(length(groups)) gr[i] <- groups[i]
          if(length(cond)) co[i] <- cond[i]
          at[[w]]     <- if(ll) 1:ll else pretty(z[i])
          labels[[w]] <- if(ll) perhapsAbbrev(l) else format(at[[w]])
          limits[[w]] <- if(ll) c(2/3, ll+1/3) else range(z[i])
        }
      if(length(cond))
        {
          nuc <- length(levels(cond))
          at <- at[rep(seq(1, length(at)), each=nuc)]
          labels <- labels[rep(seq(1, length(labels)), each=nuc)]
          limits <- limits[rep(seq(1, length(limits)), each=nuc)]
          levs <- levs[rep(seq(1, length(levs)), each=nuc)]
          
          formula <- if(!conf.int) x$yhat ~ xp | cond*p else
          Cbind(x$yhat, x$lower, x$upper) ~ xp | cond*p
        }
      else
        {
          formula <- if(!conf.int) x$yhat ~ xp | p else
          Cbind(x$yhat, x$lower, x$upper) ~ xp | p
        }
      
      pan <- function(x, y, ...)
        {
          pn <- panel.number()
          lev    <- levs[[pn]]
          col <- trellis.par.get('superpose.line')$col
          if(!length(lev) && length(unique(x[!is.na(x)])) > nlevels)
            {
              yy <- y
              if(length(perim))
                {
                  j <- perim(x, y)
                  yy[j] <- NA
                  if(length(attr(yy,'other')))
                    attr(yy, 'other')[j,] <- NA
                }
              panel.xYplot(x, yy, ...)
              if(length(data) && length(xd <- data[[names(levs)[pn]]]))
                {
                  xd <- xd[!is.na(xd)]
                  doscat1d(xd, approx(x, y, xout=xd, rule=2, ties=mean)$y,
                         col=col[1])
                }
            }
          else
            {
              panel.points(x, y, pch=19)
              yoth <- attr(y, 'other')
              if(length(yoth)) for(u in unique(x))
                llines(c(u,u), yoth[x==u,])
            }
          addpanel(x, y, ...)
        }
      scales <- list(x=list(relation='free', limits=limits,
                       at=at, labels=labels))
      r <- list(formula=formula, groups=gr, subset=subset, type='l',
                method=if(conf.int) 'filled bands' else 'bars',
                col.fill=col.fill,
                xlab='', ylab=ylab, ylim=ylim,
                panel=pan, scales=scales)
      if(length(dotlist)) r <- c(r, dotlist)
      if(length(sub   )) r$sub    <- sub
    }
  else
    {
      v <- character(0)
      bar <- ''
      f <- if(!missing(formula)) gsub(' .*','',as.character(formula)[2])
      else varying[1]
      iv <- var.inner(as.formula(paste('~', f)))
      if(missing(xlab)) xlab <- labelPlotmath(label[iv], units[iv])
      if(missing(formula))
        {
          xvar <- varying[1]
          ## change formula like ~x|foo to x
          if(length(varying) > 1)
            {
              v <- varying[-1]
              if(length(gname))
                {
                  groups <- x[[gname]]
                  v <- setdiff(v, gname)
                }
              else
                {
                  nu <- sapply(x[v], function(u) length(unique(u)))
                  if(!predpres && any(nu <= nlevels))
                    {
                      i <- which.min(nu)
                      gname <- v[i]
                      groups <- x[[gname]]
                      v <- setdiff(v, gname)
                    }
                }
              if(length(v))
                {
                  bar <- paste(v, collapse='*')
                  for(z in v)
                    {
                      if(z %in% names(x) && is.numeric(x[[z]]))
                        {
                          x[[z]] <- factor(format(x[[z]]))
                          levels(x[[z]]) <- paste(z,':',levels(x[[z]]), sep='')
                        }
                    }
                }
            } 
          xv <- x[[xvar]]
          xdiscrete <- is.factor(xv) || is.character(xv) ||
                       length(unique(xv[!is.na(xv)])) <= nlevels
          if(xdiscrete)
            {
              f <- paste(xvar, if(conf.int) 'Cbind(yhat,lower,upper)'
              else 'yhat', sep='~')
              if(bar != '') f <- paste(f, bar, sep='|')
              formula <- eval(parse(text=f))
              if(length(v)) for(z in v)
                {
                  if(z %in% names(x) && is.numeric(x[[z]]))
                    {
                      x[[z]] <- factor(format(x[[z]]))
                      levels(x[[z]]) <- paste(z,':',levels(x[[z]]), sep='')
                    }
                }
              r <- Dotplot(formula, groups=groups, subset=subset,
                           xlim=ylim, xlab=ylab, ylab=xlab,
                           sub=sub, data=x, ...)
              return(r)
            }
          if(bar != '') f <- paste(f, '|', bar)
        }
      else
        {
          f <- as.character(formula)[2]
          xvar <- gsub(' .*', '', f)
          if(length(grep('\\|', f)))
            {
              g <- gsub(' ', '', f)
              g <- gsub('.*\\|', '', g)
              v <- strsplit(g, '\\*')[[1]]
            }
          if(length(v)) for(z in v)
            {
              if(z %in% names(x) && is.numeric(x[[z]]))
                {
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
      if(xdiscrete)
        {
          xv <- as.factor(xv)
          xlev <- levels(xv)
          xscale <- list(x=list(at=1:length(xlev), labels=xlev))
          x[[xvar]] <- as.integer(xv)
        }
      
      pan <- function(x, y, groups=NULL, subscripts, ...)
        {
          yy <- y
          if(length(perim))
            {
              j <- !perim(x, y)
              yy[j] <- NA
              if(length(attr(yy,'other')))
                attr(yy, 'other')[j,] <- NA
            }
          panel.xYplot(x, yy, groups=groups, subscripts=subscripts, ...)
          col <- trellis.par.get('superpose.line')$col

          if(length(groups) && length(gd <- data[[gname]]) &&
             length(xd <- data[[xvar]]))
            {
              g <- groups[subscripts]
              j <- 0
              for(w in levels(g))
                {
                  j <- j + 1
                  z  <- g==w
                  xg <- x[z]
                  yg <- y[z]
                  x1 <- xd[gd==w]
                  x1 <- x1[!is.na(x1)]
                  doscat1d(x1, approx(xg, yg, xout=x1, rule=2, ties=mean)$y,
                           col=col[j])
                }
            }
          else if(length(xd <- data[[xvar]]))
            {
              xd <- xd[!is.na(xd)]
              doscat1d(xd, approx(x, y, xout=xd, rule=2, ties=mean)$y,
                       col=col[1])
            }
          addpanel(x, y, groups=NULL, subscripts=subscripts, ...)
        }

      r <- list(formula=formula, data=x, subset=subset,
                type=if(xdiscrete) 'b' else 'l',
                method=if(conf.int) 'filled bands' else 'bars',
                col.fill=col.fill,
                xlab=xlab, ylab=ylab, ylim=ylim, panel=pan)
      if(length(xscale)) r$scales <- xscale
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
  if(lattice)
    {
      z <- 
        function(x, y, ..., xx, yy, text, cex, adj, family)
          ltext(xx, yy, text, cex=cex, adj=adj, fontfamily=family)
      formals(z) <-
        eval(substitute(alist(x=, y=, ...=, xx=xx, yy=yy, text=k,
                              cex=cex, adj=adj, family=fam),
                        list(xx=x, yy=y, k=k, cex=cex, adj=adj,
                             fam=fam)))
    }
  else
    {
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
