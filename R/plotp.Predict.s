plotp.Predict <-
  function(data,
           subset, xlim, ylim, xlab, ylab,
           rdata=NULL, nlevels=3,
           vnames=c('labels', 'names'),
           histSpike.opts=list(frac=function(f) 0.01 + 
                                 0.02 * sqrt(f - 1)/sqrt(max(f, 2) - 1),
                               side=1, nint=100),
           ncols=3, width=800, ...)
{
  varypred <- ('.set.'       %in%  names(data)) &&
              ('.predictor.' %nin% names(data))
  if(varypred) {
    data$.predictor. <- data$.set.
    data$.set. <- NULL
  }
  predpres <- length(data$.predictor.) > 0
  vnames   <- match.arg(vnames)
  
  dohist <- function(...) {
    so <- histSpike.opts
    do.call('histSpikeg', c(list(...), so))
  }

  info     <- attr(data, 'info')
  at       <- info$Design
  label    <- at$label
  units    <- at$units
  adjust   <- info$adjust
  varying  <- setdiff(info$varying, '.set.')
  if(predpres && identical(sort(unique(data$.predictor.)), sort(varying)))
    varying <- NULL
  conf.int <- info$conf.int

  if(length(varying) > 2)
    stop('more than 2 varying variables not allowed')

  pmlabel <- character(length(label))
  names(pmlabel) <- names(label)
  for(i in 1 : length(label))
    pmlabel[i] <-markupSpecs$html$varlabel(label[i], units[i])
  
  if(predpres)
    data$.Predictor. <-
      switch(vnames,
             names  = data$.predictor.,
             labels = pmlabel[as.character(data$.predictor.)] )
    
  if(! missing(subset)) {
    subset <- eval(substitute(subset), data)
    data   <- data[subset,, drop=FALSE]
  }

  if(missing(ylab)) ylab <- info$ylabhtml
  
  if(! length(data$lower)) conf.int <- FALSE

  cllab <- if(conf.int) paste(conf.int, 'C.L.')
  
  if(missing(ylim))
    ylim <- if(conf.int) 
              with(data, c(min(c(yhat, lower), na.rm=TRUE),
                           max(c(yhat, upper), na.rm=TRUE)))
            else
              range(data$yhat, na.rm=TRUE)

  adjto <- paste0('Adjusted to:<br>', adjust)
  if(predpres) names(adjto) <- unique(data$.predictor.)

  fm <- function(x) format(x, digits=4)
  
  if(predpres) {   ## User did not specify which predictors to plot; all plotted
    data$.predictor.  <- factor(data$.predictor., unique(data$.predictor.))

    ## Determine which predictors are discrete
    isdiscrete <- function(z) is.factor(z) || is.character(z) ||
        length(unique(z[!is.na(z)])) <= nlevels
    lp    <- levels(data$.predictor.)
    isdis <- sapply(data[lp], isdiscrete)

    ## Do all continuous predictors
    vcon  <- lp[! isdis]
    ncont <- 0
    cont  <- list()
    height <- 400 * ceiling(length(vcon) / ncols)

    for(v in vcon) {
      ncont <- ncont + 1
      dat <- data[data$.predictor. == v,, drop=FALSE]
      dat$.x. <- dat[[v]]
      xlab <- pmlabel[v]
      ht   <- with(dat, paste0(v, '=', fm(.x.), '<br>',
                               fm(yhat), ' [', fm(lower), ',',
                               fm(upper), ']'))

      if(length(varying) != 2) {
        ht[1] <- paste0(ht[1], '<br>', adjto[v])
        dat$.ht. <- ht
        a <- plotly::plot_ly(dat, height=height, width=width)
        a <- plotly::add_lines(a, x=~.x., y=~yhat, text=~.ht., color=I('black'),
                               hoverinfo='text',
                               name='Estimate', legendgroup='Estimate',
                               showlegend=ncont == 1)
        if(conf.int)
          a <- plotly::add_ribbons(a, x=~.x., ymin=~lower, ymax=~upper,
                                   color=I('lightgray'), hoverinfo='none',
                                   name=cllab, legendgroup=cllab,
                                   showlegend=ncont == 1)
        if(length(rdata) && v %in% names(rdata)) {
          form <- as.formula(paste('yhat ~', v))
          a <- histSpikeg(form, data=rdata, predictions=dat, ylim=ylim,
                          plotly=a, showlegend=ncont == 1)
        }
      } else { # a second variable (for superpositioning) is varying
        w <- varying[2]
        dat$.g. <- dat[[w]]
        j <- which(dat$.x. == min(dat$.x.))
        ht[j] <- paste0(ht[j], '<br>', adjto[v])
        dat$.ht. <- ht
        a <- plotly::plot_ly(dat, height=height, width=width)
        a <- plotly::add_lines(a, x=~.x., y=~yhat, text=~.ht., color=~.g.,
                               hoverinfo='text',
                               name='Estimate', legendgroup='Estimate',
                               showlegend=ncont == 1)
        if(conf.int)
          a <- plotly::add_ribbons(a, x=~.x., ymin=~lower, ymax=~upper,
                                   color=~.g., hoverinfo='none',
                                   name=cllab, legendgroup=cllab,
                                   showlegend=ncont == 1)
        if(length(rdata) && all(c(v, w) %in% names(rdata))) {
          form <- as.formula(paste('yhat ~', v, '+', w))
          a <- histSpikeg(form, data=rdata, predictions=dat, ylim=ylim,
                          plotly=a, showlegend=ncont == 1)
        }
      }
      a <- plotly::layout(a,
                          xaxis=list(title=xlab),
                          yaxis=list(title=ylab))
      cont[[ncont]] <- a
    }
    if(ncont > 0) {
      if(ncont == 1) cont <- cont[[1]]
      else {
        nrows <- ceiling(ncont / ncols)
        cont <- plotly::subplot(cont, nrows=nrows, shareY=TRUE, titleX=TRUE)
        }
    }
    
    ## Do all categorical predictors
    if(sum(isdis) == 0) return(cont)
    
    vcat  <- lp[isdis]
    ncat  <- 0
    catg  <- list()
    nlev  <- integer(length(vcat))

    major <- minor <- character(0)
    X <- Lower <- Upper <- numeric(0)
    for(v in vcat) {
      ncat <- ncat + 1
      dat <- data[data$.predictor. == v,, drop=FALSE]
      dat$.x. <- dat[[v]]
      xlab <- pmlabel[v]

      X <- c(X, dat$yhat)
      if(conf.int) {
        Lower <- c(Lower, dat$lower)
        Upper <- c(Upper, dat$upper)
      }
      minor <- c(minor, as.character(dat[[v]]))
      major <- c(major, rep(xlab, nrow(dat)))
    }

    catg <- dotchartpl(X, major, minor, lower=Lower, upper=Upper,
                       htext=format(X, digits=4), xlab=ylab,
                       tracename='Estimate', limitstracename=cllab,
                       width=width)
    return(list(Continuous=cont, Categorical=catg))
  }

  ## .predictor. not present; assume one plot

  v        <- varying[1]
  data$.x. <- data[[v]]
  if(length(varying) > 1) {
    w <- varying[2]
    data$.g. <- data[[w]]
    }
  ht <- with(data, paste0(v, '=', fm(data$.x.), '<br>',
                          fm(yhat), ' [', fm(lower), ', ',
                          fm(upper), ']'))
  j <- which(data$.x. == min(data$.x.))
  ht[j] <- paste0(ht[j], '<br>', adjto)
  data$.ht. <- ht

  a <- plotly::plot_ly(data)
  if(length(varying) == 1) {
    a <- plotly::add_lines(a, x=~.x., y=~yhat, color=I('black'),
                           text=~.ht., hoverinfo='text',
                           name='Estimate')
    a <- plotly::add_ribbons(a, x=~.x., ymin=~lower, ymax=~upper,
                             hoverinfo='none', name=cllab,
                             color=I('lightgray'))
    if(length(rdata) && varying %in% names(rdata)) {
      form <- as.formula(paste('yhat ~', v))
      a <- histSpikeg(form, predictions=data, data=rdata,
                      plotly=a, ylim=ylim)
    }
  } else {  # superpositioning (grouping) variable also present
    a <- plotly::add_lines(a, x=~.x., y=~yhat, color=~.g.,
                           text=~.ht., hoverinfo='text')
    a <- plotly::add_ribbons(a, x=~.x., ymin=~lower, ymax=~upper,
                             color=~.g., hoverinfo='none')
    if(length(rdata) && all(varying %in% names(rdata))) {
      form <- as.formula(paste('yhat ~', v, '+', w))
      a <- histSpikeg(form, predictions=data, data=rdata,
                      plotly=a, ylim=ylim)
      }
  }

  if(missing(xlab)) xlab <- pmlabel[v]
  if(missing(xlim)) xlim <- NULL  #range(data$.x.)
  plotly::layout(a, xaxis=list(title=xlab, range=xlim),
                    yaxis=list(title=ylab, range=ylim))
}
