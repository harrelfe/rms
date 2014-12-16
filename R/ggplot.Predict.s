ggplot.Predict <-
  function(data, formula, groups=NULL,
           aestype=c('color', 'linetype'),
           varypred=FALSE, subset,
           xlim, ylim, xlab, ylab,
           rdata=NULL, anova=NULL, pval=FALSE, size.anova=4,
           adj.subtitle, size.adj=2.5, perim=NULL, nlevels=3,
           flipxdiscrete=TRUE,
           legend.position='right', layout=NULL, addlayer=NULL,
           histSpike.opts=list(frac=0.02, side=1, nint=100), type=NULL,
           ...)
{
  if(varypred) {
    data$.predictor. <- data$.set.
    data$.set. <- NULL
  }
  predpres <- length(data$.predictor.) > 0

  if(predpres && missing(legend.position)) legend.position <- 'top'
  
  dohist <- function(...) {
    so <- histSpike.opts
    ##    if(!length(so$col)) so$col <- col
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
  dotlist  <- list(...)

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
    tanova <- function(name, x, y)
      annotateAnova(name, stat, x, y, ggplot=TRUE, xlim=xlim, ylim=ylim,
                         size=size.anova)
  } else tanova <- function(...) {}

  ## See http://bigdata-analyst.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
  ## Changed to expect size in mm
  footnote <- function(text, size=2.5, color=grey(.5)) {
    grid::pushViewport(viewport())
    grid::grid.text(label= text,
                    x   = unit(1,"npc") - unit(2, "mm"),
                    y   = unit(2, "mm"),
                    just= c("right", "bottom"),
                    gp  = gpar(fontsize=size/0.3527778, col=color))
   grid::popViewport()
}
  if(predpres) {
    if(! missing(formula))
      stop('formula may not be given when multiple predictors were varied separately')
    grid.newpage()
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
    pushViewport(viewport(layout = grid.layout(layout[1], layout[2])))

    .co <- if(length(groups))  as.factor(data[[groups]])
    nr <- 1; nc <- 0
    for(w in lp) {
      nc <- nc + 1
      if(nc > layout[2]) {nr <- nr + 1; nc <- 1}
      i  <- p == w
      z    <- data[i, w]
      l  <- levels(z)
      ll <- length(l)
      xlim <- if(ll) c(1, ll) else range(pretty(z))
      yhat <- data[i, 'yhat']
      xl <- labelPlotmath(label[w], units[w])
      zz <- data.frame(.z=z, .yhat=yhat)
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
      if(length(type))
        g <- g + switch(type, p=geom_point(), l=geom_line(),
                        b=geom_point() + geom_line())
      else 
        g <- g + if(ll) geom_point() else geom_line()
      ## For unknown reasons + xlab(xl) yields "argument xlab missing"
      g <- g + ylim(ylim) + labs(x=xl, y=ylab) +
        theme(plot.margin = grid::unit(rep(.2, 4), 'cm'))
      if(length(groups)) {
        if(nr == 1 && nc == 1) {
          f <- get(paste('scale', aestype[1], 'discrete', sep='_'))
          labl <- labelPlotmath(label[groups], units[groups])
          g <- g + if(aestype[1] == 'size') f(name=labl, range=c(.2, 1.5)) else
                   f(name=labl)
          g <- g + theme(legend.position=legend.position)
        } else g <- g + theme(legend.position='none')
      }
      g <- g + if(conf.int) with(zz, tanova(w, c(.z, .z, .z),
                                   c(.yhat, lower, upper)))
       else with(zz, tanova(w, .z, yhat))
      if(length(addlayer)) {
        ## Can't specify addlayer = geom_x() + geom_x():
        ## non-numeric argument to binary operator
        if(is.list(addlayer))
          for(j in 1 : length(addlayer)) g <- g + addlayer[[j]]
        else g <- g + addlayer
      }
      if(conf.int)
        g <- g +
        if(ll) geom_errorbar(aes(ymin=lower, ymax=upper), width=0)
        else geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2,
                         linetype=0, show_guide=FALSE)
      if(length(rdata) && w %in% names(rdata)) {
        rdata$.z <- rdata[[w]]
        if(length(.co)) {
          rdata$.cond <- rdata[[groups]]
          form <- .yhat ~ .z + .cond
        } else form <- .yhat ~ .z
        g <- g + dohist(form, predictions=zz, data=rdata, ylim=ylim)
      }
      print(g, vp = viewport(layout.pos.row=nr, layout.pos.col=nc))
    }
    if(length(sub)) footnote(sub, size=size.adj)
    return(invisible())

  } else  { # .predictor. not included; user specified predictors to show

    v  <- varying
    xn <- v[1]    ## name of x-axis variable (first variable given to Predict)
    if(missing(xlab)) xlab <- labelPlotmath(label[xn], units[xn])
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
    if(xdiscrete) {
      if(flipxdiscrete) g <- g + coord_flip()
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
          f <- get(paste('scale', aestype[j], 'discrete', sep='_'))
          labl <- labelPlotmath(label[groups[j]], units[groups[j]])
          g <- g + if(aestype[j] == 'size') f(name=labl, range=c(.2, 1.5)) else
                   f(name=labl)
        }
        g <- g + theme(legend.position=legend.position)
      }
      if(conf.int) g <- g + geom_ribbon(data=data, aes(ymin=lower, ymax=upper),
                                        alpha=0.2, linetype=0,
                                        show_guide=FALSE)
    }
    g <- g + if(conf.int) tanova(xn, c(xv, xv, xv),
                                 c(data$yhat, data$lower, data$upper)) else
     tanova(xn, xv, data$yhat)
                    
    if(length(rdata) && xn %in% names(rdata)) {
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
      if(is.list(addlayer))
        for(j in 1 : length(addlayer)) g <- g + addlayer[[j]]
      else g <- g + addlayer
    }
    if(! length(sub)) return(g)
    print(g)  # needed for footnote to work
    footnote(sub, size=size.adj)
    invisible()
  }
}

utils::globalVariables(c('.z', '.yhat', 'lower', 'upper'))
