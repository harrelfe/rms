ggplot.Predict <-
  function(data, formula, groups=NULL,
           aestype=c('color', 'linetype'),
           conf=c('fill', 'lines'),
           varypred=FALSE, subset,
           xlim, ylim, xlab, ylab,
           colorscale=function(...)
             scale_color_manual(..., values=c("#000000", "#E69F00", "#56B4E9",
               "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")),
           rdata=NULL, anova=NULL, pval=FALSE, size.anova=4,
           adj.subtitle, size.adj=2.5, perim=NULL, nlevels=3,
           flipxdiscrete=TRUE,
           legend.position='right', legend.label=NULL,
           vnames=c('labels', 'names'), abbrev=FALSE, minlength=6,
           layout=NULL, addlayer=NULL,
           histSpike.opts=list(frac=0.02, side=1, nint=100),
           type=NULL, ...)
{
  if(varypred) {
    data$.predictor. <- data$.set.
    data$.set. <- NULL
  }
  predpres <- length(data$.predictor.) > 0

  if(predpres && missing(legend.position)) legend.position <- 'top'
  conf <- match.arg(conf)
  vnames <- match.arg(vnames)
  
  dohist <- function(...) {
    so <- histSpike.opts
    do.call('histSpikeg', c(list(...), so))
  }

  
  info     <- attr(data, 'info')
  at       <- info$Design
  label    <- at$label
  units    <- at$units
  values   <- info$values
  adjust   <- info$adjust
  yunits   <- info$yunits
  varying  <- info$varying
  conf.int <- info$conf.int

  pmlabel <- vector('expression', length(label))
  names(pmlabel) <- names(label)
  for(i in 1 : length(label))
    pmlabel[i] <- labelPlotmath(label[i], units[i])

  glabel <- function(gname, j=1) {
    if(! length(legend.label)) pmlabel[gname]
    else if(is.logical(legend.label)) ''
    else legend.label[j]
  }
  
  if(! missing(subset)) {
    subset <- eval(substitute(subset), data)
    data <- data[subset,, drop=FALSE]
  }

  if(length(groups) && is.logical(groups) && ! groups) groups <- NULL
  else if(length(groups)) {
    if(length(groups) > 2 || !is.character(groups) ||
       any(groups %nin% names(data)))
      stop('groups must be one or two predictor names')
  } else if(! predpres && length(varying) > 1) groups <- varying[2]

  ## Make all grouping variables discrete for proper aesthetic mapping
  if(length(groups)) for(v in groups) data[[v]] <- as.factor(data[[v]])
  
  if(missing(ylab))     ylab     <- info$ylabPlotmath
  if(! length(data$lower)) conf.int <- FALSE
  
  if(missing(ylim))
    ylim <- range(pretty(
             if(conf.int) c(data$yhat, data$lower, data$upper)
              else data$yhat), na.rm=TRUE)

  if(missing(adj.subtitle)) adj.subtitle <- length(adjust) > 0
  sub <- if(adj.subtitle && length(adjust)==1)
    paste('Adjusted to:', adjust, sep='') else NULL
  
  if(length(anova)) {
    stat   <- plotmathAnova(anova, pval)
    tanova <- function(name, x, y, flip=FALSE)
      if(flip)
        annotateAnova(name, stat, y, x, ggplot=TRUE, xlim=ylim, ylim=xlim,
                      size=size.anova) else
        annotateAnova(name, stat, x, y, ggplot=TRUE, xlim=xlim, ylim=ylim,
                      size=size.anova)
  } else tanova <- function(...) {}

  ## See http://bigdata-analyst.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
  ## size is in mm
  footnote <- function(object, text, size=2.5, color=grey(.5))
    arrangeGrob(object, sub = textGrob(text, x = 1, hjust = 1.01,
       vjust=0.1, gp = gpar(fontsize =size/0.3527778 )))
  
  if(predpres) {
    p  <- as.factor(data$.predictor.)
    xp <- rep(NA, length(p))
    levs <- at <- labels <- limits <- list()
    lp <- levels(p)
    np <- length(lp)
    if(! length(layout)) 
      layout <-
           if(np <= 4) c(2,2)
      else if(np <= 6) c(2,3)
      else if(np <= 9) c(3,3)
      else if(np <=12) c(3,4)
      else if(np <=16) c(4,4)
      else             c(4,5)
#    pushViewport(viewport(layout = grid.layout(layout[1], layout[2])))

    .co <- if(length(groups))  as.factor(data[[groups]])
    # nr <- 1; nc <- 0
    Plt <- list()
    for(w in lp) {
      # nc <- nc + 1
      # if(nc > layout[2]) {nr <- nr + 1; nc <- 1}
      i  <- p == w
      z    <- data[i, w]
      l  <- levels(z)
      if(abbrev) {
        l <- abbreviate(l, minlength=minlength)
        levels(z) <- l
      }
      ll <- length(l)
      xlim <- if(ll) c(1, ll) else range(pretty(z))
      yhat <- data[i, 'yhat']
      xl <- if(vnames == 'labels') pmlabel[w] else w
      zz <- data.frame(.z=z, .yhat=yhat)
      if(! missing(formula))
        zz <- cbind(zz, data[i, all.vars(formula), drop=FALSE])
      if(conf.int) {
        zz$lower <- data[i, 'lower']
        zz$upper <- data[i, 'upper']
      }
      if(length(.co)) {
        zz$.cond <- .co[i]
        g <- eval(parse(text=sprintf(
           'ggplot(zz, aes(x=.z, y=.yhat, %s=.cond))', aestype[1])))
      }
      else g <- ggplot(zz, aes(x=.z, y=.yhat))

      xdiscrete <- is.factor(z) || is.character(z) ||
                   length(unique(z[!is.na(z)])) <= nlevels
      flipped <- FALSE
      if(xdiscrete) {
        if(flipxdiscrete && ! is.numeric(z)) {
          g <- g + coord_flip()
          flipped <- TRUE
        }
        g <- g + geom_point() +
          if(length(type) && type %in% c('l', 'b')) geom_line()
        if(is.numeric(z)) g <- g + scale_x_discrete(breaks=unique(z))
      } else {
        if(length(type))
          g <- g + switch(type, p=geom_point(), l=geom_line(),
                          b=geom_point() + geom_line())
        else g <- g + geom_line()
      }

      ## Need to following or geom_ribbon will improperly clip regions
      if(flipped) g <- g + ylim(ylim) else
                  g <- g + coord_cartesian(ylim=ylim)
      g <- g + labs(x=xl, y=ylab) +
        theme(plot.margin = grid::unit(rep(0, 4), 'cm'))
      ## use rep(.1, 4) if using print(..., viewport=...) for multiple plots
      if(length(groups)) {
        # if(nr == 1 && nc == 1) {
        if(w == lp[1]) {
          f <- if(aestype[1] == 'color') colorscale else
           get(paste('scale', aestype[1], 'discrete', sep='_'))
          labl <- glabel(groups[1])
          g <- g + if(aestype[1] == 'size') f(name=labl, range=c(.2, 1.5)) else
                   f(name=labl)
          g <- g + theme(legend.position=legend.position)
        } else g <- g + theme(legend.position='none')
      }
      xa <- if(conf.int) c(zz$.z, zz$.z, zz$.z) else zz$.z
      ya <- if(conf.int) c(zz$.yhat, zz$lower, zz$upper) else zz$.yhat
      g <- g + tanova(w, xa, ya, flip=flipped)
      
      if(length(addlayer)) {
        ## Can't specify addlayer = geom_x() + geom_x():
        ## non-numeric argument to binary operator
        if(! length(class(addlayer)) && is.list(addlayer))
          for(j in 1 : length(addlayer)) g <- g + addlayer[[j]]
        else g <- g + addlayer
      }
      if(conf.int) {
        if(ll || xdiscrete) g <- g +
          geom_errorbar(aes(ymin=lower, ymax=upper), width=0)
        else {
          if(conf == 'fill')
            g <- g + geom_ribbon(aes(ymin=lower, ymax=upper),
                                 alpha=0.2, linetype=0,
                                 show_guide=FALSE)
          else {
            ## geom_ribbon with fill=NA draws vertical lines at
            ## ends of confidence regions
            g <- g + geom_line(aes(x=.z, y=lower))
            g <- g + geom_line(aes(x=.z, y=upper))
          }
        }
      }

      if(! missing(formula)) g <- g + facet_grid(formula)
      
      if(! is.factor(z) && length(rdata) && w %in% names(rdata)) {
        rdata$.z <- rdata[[w]]
        if(length(.co)) {
          rdata$.cond <- rdata[[groups]]
          form <- .yhat ~ .z + .cond
        } else form <- .yhat ~ .z
        g <- g + dohist(form, predictions=zz, data=rdata, ylim=ylim)
      }
      # print(g, vp = viewport(layout.pos.row=nr, layout.pos.col=nc))
      Plt[[w]] <- g
    }
    Plt <- do.call(arrangeGrob, c(Plt, list(ncol=layout[2])))
    if(length(sub)) Plt <- footnote(Plt, size=size.adj)
    return(Plt)

  } else  { # .predictor. not included; user specified predictors to show

    v  <- varying
    xn <- v[1]    ## name of x-axis variable (first variable given to Predict)
    if(missing(xlab)) xlab <- pmlabel[xn]
    xv <- data[[xn]]
    xdiscrete <- is.factor(xv) || is.character(xv) ||
                 length(unique(xv[!is.na(xv)])) <= nlevels
    if(length(perim)) {
      j <- if(! length(groups)) perim(xv, NULL) else perim(xv, data[[groups[1]]])
      data$yhat[! j] <- NA
      if(conf.int) data$lower[! j] <- data$upper[! j] <- NA
    }
    ae <- paste('aes(x=', xn, ', y=yhat', sep='')
    if(length(groups)) for(j in 1 : length(groups))
      ae <- paste(ae, ', ', aestype[j], '=', groups[j], sep='')
    ae <- eval(parse(text=paste(ae, ')', sep='')))
    class(data) <- setdiff(class(data), 'Predict')  # so won't involve ggplot.Predict
    g <- ggplot(data, ae) + labs(x=xlab, y=ylab)
    flipped <- FALSE
    if(xdiscrete) {
      if(flipxdiscrete && ! is.numeric(xv)) {
        g <- g + coord_flip()
        flipped <- TRUE
      }
      g <- g + geom_point() +
        if(length(type) && type %in% c('l', 'b')) geom_line()
      if(conf.int) g <- g + geom_errorbar(aes(ymin=lower, ymax=upper), width=0)
    } else {
      if(length(type))
        g <- g + switch(type, p=geom_point(), l=geom_line(),
                        b=geom_point() + geom_line())
      else g <- g + geom_line()
      if(length(groups)) {
        for(j in 1 : length(groups)) {
          f <- if(aestype[j] == 'color') colorscale else
               get(paste('scale', aestype[j], 'discrete', sep='_'))
          labl <- glabel(groups[j], j)
          g <- g + if(aestype[j] == 'size') f(name=labl, range=c(.2, 1.5)) else
                   f(name=labl)
        }
        g <- g + theme(legend.position=legend.position)
      }
      if(conf.int) {
        if(conf == 'fill')
          g <- g + geom_ribbon(data=data, aes(ymin=lower, ymax=upper),
                 alpha=0.2, linetype=0,
                 show_guide=FALSE)
        else {
          g <- g + eval(parse(text=
               sprintf('geom_line(data=data, aes(x=%s, y=lower))', xn)))
          g <- g + eval(parse(text=
               sprintf('geom_line(data=data, aes(x=%s, y=upper))', xn)))
        }
      }
    }
    if(flipped) g <- g + ylim(ylim) else
                g <- g +  coord_cartesian(ylim=ylim)
    xa <- if(conf.int) c(xv, xv, xv) else xv
    ya <- if(conf.int) c(data$yhat, data$lower, data$upper) else data$yhat
    g <- g + tanova(xn, xa, ya, flip=flipped)
                    
    if(! is.factor(xv) && length(rdata) && xn %in% names(rdata)) {
      form <- paste('yhat', xn, sep='~')
      if(length(groups)) form <- paste(form, groups[1], sep='+')
      form <- eval(parse(text=form))
      g <- g + dohist(form, predictions=data, data=rdata, ylim=ylim)
    }
    ## Get list of varying variables that are not a groups variable
    ## These will be for faceting
    ## If the faceting formula is specified, just use it

    f <- if(length(v) > 1) setdiff(v[-1], groups)
    if(length(f)) {
      if(missing(formula)) {
        k <- length(f)
        formula <- if(k == 1) paste('~', f[1])
        else if(k == 2) paste(f[1], f[2], sep='~')
        else if(k == 3) paste(f[1], '~', f[2], '*', f[3])
        else if(k == 4) paste(f[1], '*', f[2], '~', f[3], '*', f[4])
        else stop('too many varying variables to use faceting')
        formula <- as.formula(formula)
      }
      g <- g + facet_grid(formula)
    }
    if(length(addlayer)) {
      if(! length(class(addlayer)) && is.list(addlayer))
        for(j in 1 : length(addlayer)) g <- g + addlayer[[j]]
      else g <- g + addlayer
    }
    if(length(sub)) g <- footnote(g, sub)
    g
  }
}

utils::globalVariables(c('.z', '.yhat', 'lower', 'upper'))
