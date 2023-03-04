bplot <-
  function(x, formula, lfun=lattice::levelplot, xlab, ylab, zlab,
           adj.subtitle=!info$ref.zero, cex.adj=.75, cex.lab=1,
           perim, showperim=FALSE,
           zlim=range(yhat, na.rm=TRUE),
           scales=list(arrows=FALSE), xlabrot, ylabrot, zlabrot=90, ...)
{
  if(! requireNamespace('lattice', quietly=TRUE))
    stop('lattice package not installed')
  
  lfunname <- deparse(substitute(lfun))
  if(missing(xlabrot))
    xlabrot <- switch(lfunname,
                      wireframe=30, contourplot=0, levelplot=0, 0)
  if(missing(ylabrot))
    ylabrot <- switch(lfunname,
                      wireframe=-40, contourplot=90, levelplot=90, 0)
  info    <- attr(x, 'info')
  varying <- info$varying
  if(length(varying) < 2) stop('should vary at least two variables')
  if(missing(formula))
    {
      nx      <- varying[1]
      ny      <- varying[2]
      formula <- paste('yhat ~', nx, '*', ny)
      if(length(varying) > 2)
        formula <- paste(formula, '|', paste(varying[-(1:2)], collapse='*'))
      formula <- as.formula(formula)
    }
  else
    {
      ter <- attributes(terms(formula))
      vars <- ter$term.labels
      nx <- vars[1]
      ny <- vars[2]
      if(!ter$response)
        formula <- as.formula(paste('yhat', format(formula)))
    }

  data    <- x
  yhat    <- x$yhat
  y       <- x[[ny]]
  x       <- x[[nx]]

  at      <- info$Design
  label   <- at$label
  units   <- at$units

  if(missing(xlab)) xlab <- labelPlotmath(label[nx], units[nx])
  xlab <- list(label=xlab, rot=xlabrot, cex=cex.lab)
  if(missing(ylab)) ylab <- labelPlotmath(label[ny], units[ny])
  ylab <- list(label=ylab, rot=ylabrot, cex=cex.lab)
  if(missing(zlab)) zlab <- info$ylabPlotmath
  zlab  <- list(label=zlab, rot=zlabrot, cex=cex.lab)
  
  adjust  <- info$adjust
  
  if(!missing(perim))
    {
      Ylo <- approx(perim[,1], perim[,2], x, ties=mean)$y
      Yhi <- approx(perim[,1], perim[,3], x, ties=mean)$y
      Ylo[is.na(Ylo)] <-  1e30
      Yhi[is.na(Yhi)] <- -1e30
      yhat[y < Ylo] <- NA
      yhat[y > Yhi] <- NA
      data$yhat     <- yhat
    }
  else if(showperim) stop('cannot request showperim without specifying perim')
  
  sub <- if(adj.subtitle && length(info$adjust))
    list(label=paste('Adjusted to:', info$adjust), cex=cex.adj) else NULL

  pan <- function(...)
    {
      fname <- paste('lattice::panel', gsub('lattice::', '', lfunname),
                     sep='.')
      f <- eval(parse(text = fname))
      do.call(f, list(...))
      if(showperim)
        {
          lattice::llines(perim[,'x'], perim[,'ymin'], col=gray(.85))
          lattice::llines(perim[,'x'], perim[,'ymax'], col=gray(.85))
        }
    }
  lfun(formula, panel=pan, scales=scales, zlim=zlim, ..., data=data,
       xlab=xlab, ylab=ylab, zlab=zlab, sub=sub)
}

perimeter <- function(x, y, xinc=diff(range(x))/10, n=10,
                      lowess.=TRUE)
{

  s <- !is.na(x + y)
  x <- x[s]
  y <- y[s]
  m <- length(x)
  if(m < n)
    stop("number of non-NA x must be >= n")

  i <- order(x)
  x <- x[i]
  y <- y[i]
  s <- n:(m-n+1)
  x <- x[s]
  y <- y[s]

  x <- round(x/xinc)*xinc

  g <- function(y, n)
    {
      y <- sort(y)
      m <- length(y)
      if(n > (m - n + 1)) c(NA, NA)
      else c(y[n], y[m-n+1])
    }

  r <- unlist(tapply(y, x, g, n=n))
  i <- seq(1, length(r), by=2)
  rlo <- r[i]
  rhi <- r[-i]
  s <- !is.na(rlo + rhi)
  if(!any(s))
    stop("no intervals had sufficient y observations")

  x <- sort(unique(x))[s]
  rlo <- rlo[s]
  rhi <- rhi[s]
  if(lowess.)
    {
      rlo <- lowess(x, rlo)$y
      rhi <- lowess(x, rhi)$y
    }

  structure(cbind(x, rlo, rhi),
            dimnames=list(NULL,
              c("x","ymin","ymax")), class='perimeter')
}
