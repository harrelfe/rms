ggplot.Predict <-
  function(data, formula, groups=NULL,
           aestype=c('color', 'linetype'),
           conf=c('fill', 'lines'),
           varypred=FALSE,
           sepdiscrete=c('no', 'list', 'vertical', 'horizontal'),
           subset, xlim, ylim, xlab, ylab,
           colorscale=function(...)
             scale_color_manual(..., values=c("#000000", "#E69F00", "#56B4E9",
               "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")),
           rdata=NULL, anova=NULL, pval=FALSE, size.anova=4,
           adj.subtitle, size.adj=2.5, perim=NULL, nlevels=3,
           flipxdiscrete=TRUE,
           legend.position='right', legend.label=NULL,
           vnames=c('labels', 'names'), abbrev=FALSE, minlength=6,
           layout=NULL, addlayer=NULL,
           histSpike.opts=list(frac=function(f) 0.01 + 
                                 0.02 * sqrt(f - 1)/sqrt(max(f, 2) - 1),
             side=1, nint=100),
           type=NULL, ...)
{
  sepdiscrete <- match.arg(sepdiscrete)
  class(data) <- setdiff(class(data), 'Predict')
  ## so won't involve ggplot.Predict

  if(varypred) {
    data$.predictor. <- data$.set.
    data$.set. <- NULL
  }
  predpres <- length(data$.predictor.) > 0

  if(predpres && missing(legend.position)) legend.position <- 'top'
  conf   <- match.arg(conf)
  vnames <- match.arg(vnames)
  
  dohist <- function(...) {
    so <- histSpike.opts
    do.call('histSpikeg', c(list(...), so))
  }

  ## The following exists to nullify invisible() used in arrangeGrob's
  ## returned value
  agrob <- function(...) {
    z <- gridExtra::arrangeGrob(...)
    z
  }
  
  info     <- attr(data, 'info')
  at       <- info$Design
  label    <- at$label
  units    <- at$units
  adjust   <- info$adjust
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
    tanova <- function(name, x, y, xlim, ylim, flip=FALSE,
                       empty=FALSE, dataOnly=FALSE)
      annotateAnova(name, stat, x, y, ggplot=TRUE, xlim=xlim, ylim=ylim,
                    size=size.anova, flip=flip, empty=empty, dataOnly=dataOnly)
  } else tanova <- function(...) {}

  ## See http://bigdata-analyst.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
  ## size is in mm
  footnote <- function(object, text, size=2.5, color=grey(.5))
    agrob(object, sub = grid::textGrob(text, x = 1, hjust = 1.01,
       vjust=0.1, gp = grid::gpar(fontsize =size/0.3527778 )))
  
  if(predpres) {   ## User did not specify which predictors to plot; all plotted
    data$.predictor.  <- factor(data$.predictor.)

    if(sepdiscrete != 'no') {
      ## From http://stackoverflow.com/questions/11979017/changing-facet-label-to-math-formula-in-ggplot2
      facet_wrap_labeller <- function(gg.plot, labels=NULL) {
        ## Uses functions from gridExtra
        g <- ggplotGrob(gg.plot)
        gg <- g$grobs      
        strips <- grep("strip_t", names(gg))
        for(ii in seq_along(labels))  {
          modgrob <- grid::getGrob(gg[[strips[ii]]], "strip.text", 
                             grep=TRUE, global=TRUE)
          gg[[strips[ii]]]$children[[modgrob$name]] <-
            grid::editGrob(modgrob,label=labels[ii])
        }
        g$grobs <- gg
        class(g) = c("arrange", "ggplot", class(g)) 
        g
      }

      ## Determine which predictors are discrete
      isdiscrete <- function(z) is.factor(z) || is.character(z) ||
        length(unique(z[!is.na(z)])) <= nlevels
      lp   <- setdiff(levels(data$.predictor.), groups)
      isdis <- sapply(data[lp], isdiscrete)

      dogroup <- function(type) {
        lim <- ggplot2:::limits
        v <- if(type == 'continuous') names(isdis)[! isdis] else
             names(isdis)[isdis]
        # dat <- subset(data, .predictor. %in% v)  ## would not work
        dat <- data[data$.predictor. %in% v,, drop=TRUE]
        p <- dat$.predictor.
        xx <- switch(type, continuous= numeric(nrow(dat)),
                           discrete  = character(nrow(dat)))
        for(iv in v) {
          j <- which(p == iv)
          xx[j] <- switch(type, continuous= dat[j, iv],
                                discrete  = as.character(dat[j, iv]))
        }
        dat$.xx. <- xx
        if(length(groups)) dat$.co. <- as.factor(dat[[groups]])

        if(type == 'continuous') {
          if(length(groups)) g <- eval(parse(text=
             sprintf('ggplot(dat, aes(x=.xx., y=yhat, %s=%s))',
             aestype[1], groups[1]))) + labs(x=NULL, y=ylab) + lim(ylim, 'y')
        else
          g <- ggplot(dat, aes(x=.xx., y=yhat)) + labs(x=NULL, y=ylab) +
            lim(ylim, 'y')

          g <- g + (if(length(layout))
                     facet_wrap(~ .predictor., scales='free_x',
                                ncol=layout[2]) else
                    facet_wrap(~ .predictor., scales='free_x')) + geom_line()
          if(conf.int) {
            if(conf == 'fill')
              g <- g +
                geom_ribbon(aes(x=.xx., ymin=lower, ymax=upper),
                            alpha=0.2, linetype=0, show_guide=FALSE)
            else
              g <- g +
                geom_line(aes(x=.xx., y=lower)) +
                geom_line(aes(x=.xx., y=upper))
          }

          if(length(rdata)) {
            rv <- intersect(v, names(rdata))
            rdata <- rdata[c(rv, groups)]
            ## For each variable in rdata that is in dat, set values
            ## outside the range in dat to NA.  Otherwise x-axes will
            ## be rescaled to include all raw data values, not just
            ## points at which predictions are made
            for(vv in rv) {
              a <- dat[[vv]]
              if(is.numeric(a)) {
                r <- range(a, na.rm=TRUE)
                b <- rdata[[vv]]
                i <- b < r[1] | b > r[2]
                if(any(i)) {
                  b[i] <- NA
                  rdata[[vv]] <- b
                }
              }
            }
            
            ## Reshape rdata to be tall and thin
            rdata <- reshape(as.data.frame(rdata),
                             direction='long', v.names='.xx.',
                             timevar='.predictor.',
                             varying=rv, times=rv)
            form <- 'yhat ~ .xx. + .predictor.'
            if(length(groups))
              form <- paste(form, '+', paste(groups, collapse='+'))
            form <- as.formula(form)
            g <- g + dohist(form, predictions=dat, data=rdata, ylim=ylim)
          }
        } else {   # discrete x
          if(length(groups)) g <- eval(parse(text=
                     sprintf('ggplot(dat, aes(x=yhat, y=.xx., %s=%s))',
                     aestype[1], groups[1]))) + labs(x=ylab, y=NULL)
          else
            g <- ggplot(dat, aes(x=yhat, y=.xx.)) + labs(x=ylab, y=NULL)
          g <- g + lim(ylim, 'x') +
               facet_wrap(~ .predictor., scales='free_y') + geom_point()
          if(conf.int) g <- g +
            geom_errorbarh(aes(y=.xx., xmin=lower, xmax=upper), height=0)
        }
        ## anova annotations need to be created for all variables being
        ## plotted with faceting, and annotation information must be
        ## based on a dataset with the information and the .predictor.
        ## variable, and geom_text() must be used instead of annotate()
        ## See http://stackoverflow.com/questions/2417623/manual-annotate-a-ggplot-with-different-labels-in-different-facets
        if(length(anova)) {
          .xx. <- yhat <- .label. <- hjust <- vjust <- NULL
          for(iv in v) {
            j <- which(data$.predictor. == iv)
            dat <- data[j,, drop=FALSE]
            xv <- dat[, iv]
            xx <- switch(type, continuous = xv,
                               discrete   = as.numeric(xv))
            yy <- dat[, 'yhat']
            if(conf.int) {
              xx <- c(xx, xx, xx)
              yy <- c(yy, dat[, 'lower'], dat[, 'upper'])
            }
            xlim <- if(is.factor(xv)) c(1, length(levels(xv)))
            else range(pretty(xv))
            tan <- tanova(iv, xx, yy, xlim, ylim, dataOnly=TRUE,
                          flip=type=='discrete', empty=type == 'discrete')
            .xx. <- c(.xx., tan$x)
            yhat <- c(yhat, tan$y)
            .label. <- c(.label., tan$label)
            hjust <- c(hjust, tan$hjust)
            vjust <- c(vjust, tan$vjust)
          }
          .anova. <- data.frame(.predictor.=v, .xx., yhat, .label.,
                                hjust, vjust)
          g <- g + geom_text(aes(label=.label., hjust=hjust, vjust=vjust),
                             size=size.anova,
                             data=.anova., parse=TRUE, show_guide=FALSE)
        }
        if(length(addlayer)) {
          ## Can't specify addlayer = geom_x() + geom_x():
          ## non-numeric argument to binary operator
          if(! length(class(addlayer)) && is.list(addlayer))
            for(j in 1 : length(addlayer)) g <- g + addlayer[[j]]
          else g <- g + addlayer
        }
        
        if(vnames == 'labels') g <- facet_wrap_labeller(g, pmlabel[v])
        g
    }
      
      gcont <- if(any(! isdis)) dogroup('continuous')
      gdis  <- if(any(  isdis)) dogroup('discrete')
      r <- mean(! isdis)
      return(if(length(gcont) && length(gdis))
        switch(sepdiscrete,
               list = list(continuous=gcont, discrete=gdis),
               vertical = agrob(gcont, gdis,         heights=c(r, 1-r)),
               horizontal=agrob(gcont, gdis, nrow=1, widths =c(r, 1-r)))
      else if(length(gcont)) gcont else gdis)
    }  # end if(sepdiscrete)
    ## Form separate plots and combine at end
    p <- data$.predictor.
    levs <- at <- labels <- limits <- list()
    lp   <- setdiff(levels(p), groups)
    np   <- length(lp)
    .co  <- if(length(groups))  as.factor(data[[groups]])

    if(! length(layout)) 
      layout <-
           if(np <= 4) c(2,2)
      else if(np <= 6) c(2,3)
      else if(np <= 9) c(3,3)
      else if(np <=12) c(3,4)
      else if(np <=16) c(4,4)
      else             c(4,5)
#    pushViewport(viewport(layout = grid.layout(layout[1], layout[2])))

    Plt <- list()
    jplot <- 0
    for(w in lp) {
      jplot <- jplot + 1
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
      zz <- data.frame(.xx.=z, .yhat=yhat)
      if(! missing(formula))
        zz <- cbind(zz, data[i, all.vars(formula), drop=FALSE])
      if(conf.int) {
        zz$lower <- data[i, 'lower']
        zz$upper <- data[i, 'upper']
      }
      if(length(.co)) {
        zz$.cond <- .co[i]
        g <- eval(parse(text=sprintf(
           'ggplot(zz, aes(x=.xx., y=.yhat, %s=.cond))', aestype[1])))
      }
      else g <- ggplot(zz, aes(x=.xx., y=.yhat))

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
        if(jplot == 1) {
          f <- if(aestype[1] == 'color') colorscale else
           get(paste('scale', aestype[1], 'discrete', sep='_'))
          labl <- glabel(groups[1])
          g <- g + if(aestype[1] == 'size') f(name=labl, range=c(.2, 1.5)) else
                   f(name=labl)
          g <- g + theme(legend.position=legend.position)
        } else g <- g + theme(legend.position='none')
      }
      xa <- if(conf.int) c(zz$.xx., zz$.xx., zz$.xx.) else zz$.xx.
      ya <- if(conf.int) c(zz$.yhat, zz$lower, zz$upper) else zz$.yhat
      g  <- g + tanova(w, xa, ya, xlim, ylim, flip=FALSE)   ## was flip=flipped
      
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
            g <- g + geom_line(aes(x=.xx., y=lower))
            g <- g + geom_line(aes(x=.xx., y=upper))
          }
        }
      }

      if(! missing(formula)) g <- g + facet_grid(formula)
      
      if(! is.factor(z) && length(rdata) && w %in% names(rdata)) {
        rdata$.xx. <- rdata[[w]]
        if(length(.co)) {
          rdata$.cond <- rdata[[groups]]
          form <- .yhat ~ .xx. + .cond
        } else form <- .yhat ~ .xx.
        g <- g + dohist(form, predictions=zz, data=rdata, ylim=ylim)
      }
      # print(g, vp = viewport(layout.pos.row=nr, layout.pos.col=nc))
      Plt[[jplot]] <- g
    }
    Plt <- if(jplot == 1) Plt[[1]]
    else
      do.call(agrob, c(Plt, list(ncol=layout[2])))
    if(length(sub)) Plt <- footnote(Plt, sub, size=size.adj)
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
    if(missing(xlim))
      xlim <- if(is.factor(xv)) c(1 , length(levels(xv)))
    else
      range(pretty(xv))

    g <- g + tanova(xn, xa, ya, xlim, ylim, flip=FALSE)  # was flip=flipped
                    
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

utils::globalVariables(c('.xx.', '.yhat', 'lower', 'upper'))
