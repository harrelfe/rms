ggplot.Predict <-
  function(data, formula, groups=NULL,
           aestype=c('color', 'linetype'),
           conf=c('fill', 'lines'),
           varypred=FALSE,
           sepdiscrete=c('no', 'list', 'vertical', 'horizontal'),
           subset, xlim., ylim., xlab, ylab,
           colorscale=function(...)
             scale_color_manual(..., values=c("#000000", "#E69F00", "#56B4E9",
               "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")),
           colfill='black',
           rdata=NULL, anova=NULL, pval=FALSE, size.anova=4,
           adj.subtitle, size.adj=2.5, perim=NULL, nlevels=3,
           flipxdiscrete=TRUE,
           legend.position='right', legend.label=NULL,
           vnames=c('labels', 'names'), abbrev=FALSE, minlength=6,
           layout=NULL, addlayer,
           histSpike.opts=list(frac=function(f) 0.01 + 
                                 0.02 * sqrt(f - 1)/sqrt(max(f, 2) - 1),
             side=1, nint=100),
           type=NULL, ggexpr=FALSE, ...)
{
  # .xlim, .ylim instead of xlim, ylim to distinguish from ggplot functions
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

  maddlayer <- missing(addlayer)
  if(! maddlayer) addlayer <- paste(deparse(substitute(addlayer)), collapse=' ')

  ribbonargs <- sprintf("alpha=0.2, linetype=0, fill=I('%s'), show_guide=FALSE",
                        colfill)
  
  dohist <- function(...) {
    so <- histSpike.opts
    do.call('histSpikeg', c(list(...), so))
  }

  info     <- attr(data, 'info')
  at       <- info$Design
  label    <- at$label
  units    <- at$units
  adjust   <- info$adjust
  varying  <- info$varying
  conf.int <- info$conf.int

  pmlabel <- character(length(label))
  names(pmlabel) <- names(label)
  for(i in 1 : length(label))
    pmlabel[i] <- labelPlotmath(label[i], units[i], chexpr=TRUE)

  glabel <- function(gname, j=1) {
    if(! length(legend.label)) pmlabel[gname]
    else if(is.logical(legend.label)) ''
    else legend.label[j]
  }

  ## Function to create expression( ) or "" depending on argument
  expch <- function(x) if(! length(x)) 'NULL' else
           if(grepl('expression\\(', x)) x else deparse(x)

  ## Function to construct xlim() or ylim() call
  limc <- function(limits, which)
    sprintf("%slim(%s, %s)", which, limits[1], limits[2])
  
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
  
  if(missing(ylim.))
    ylim. <- range(pretty(
      if(conf.int) c(data$yhat, data$lower, data$upper)
      else data$yhat), na.rm=TRUE)

  if(missing(adj.subtitle)) adj.subtitle <- length(adjust) > 0
  sub <- if(adj.subtitle && length(adjust)==1)
           paste('Adjusted to:', adjust, sep='') else NULL
  
  tanova <- if(length(anova))
    function(name, x, y, xlim, ylim, flip=FALSE,
             empty=FALSE, dataOnly=FALSE)
      annotateAnova(name, plotmathAnova(anova, pval),
                    x, y, ggplot=TRUE, xlim=xlim, ylim=ylim,
                    size=size.anova, flip=flip, empty=empty, dataOnly=dataOnly)
    else function(...) {}

  ## See http://bigdata-analyst.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
  ## size is in mm
  footnote <- function(object, text, size=2.5, color=grey(.5))
    arrGrob(object, sub = grid::textGrob(text, x = 1, hjust = 1.01,
          vjust=0.1, gp = grid::gpar(fontsize =size/0.3527778 )))
  
  if(predpres) {   ## User did not specify which predictors to plot; all plotted
    data$.predictor.  <- factor(data$.predictor.)

    if(sepdiscrete != 'no') {
      ## From http://stackoverflow.com/questions/11979017/changing-facet-label-to-math-formula-in-ggplot2
      ## Changed to assume that each element of labels is a character string
      ## of the form "expression(....)"
      facet_wrap_labeller <- function(gg.plot, labels=NULL) {
        ## Uses functions from gridExtra
        g <- ggplotGrob(gg.plot)
        gg <- g$grobs      
        strips <- grep("strip_t", names(gg))
        for(ii in seq_along(labels))  {
          modgrob <- grid::getGrob(gg[[strips[ii]]], "strip.text", 
                                   grep=TRUE, global=TRUE)
          gg[[strips[ii]]]$children[[modgrob$name]] <-
            grid::editGrob(modgrob,label=eval(parse(text=labels[ii])))
        }
        g$grobs <- gg
#        class(g) = c("arrange", "ggplot", class(g))
        g
      }
      
      ## Determine which predictors are discrete
      isdiscrete <- function(z) is.factor(z) || is.character(z) ||
        length(unique(z[!is.na(z)])) <= nlevels
      lp   <- setdiff(levels(data$.predictor.), groups)
      isdis <- sapply(data[lp], isdiscrete)
      
      dogroup <- function(type) {
        v <- if(type == 'continuous') names(isdis)[! isdis] else
              names(isdis)[isdis]
        # dat <- subset(data, .predictor. %in% v)  ## would not work
        dat <- data[data$.predictor. %in% v,, drop=TRUE]
        p <- dat$.predictor.
        xx <- switch(type, continuous= numeric(nrow(dat)),
                           discrete  = character(nrow(dat)))
        ## Prepare to create a "super factor" variable by concatenating
        ## all levels of all categorical variables keeping original orders
        Lev <- character()
        for(iv in v) {
          j <- which(p == iv)
          datj <- dat[j, iv]
          if(type == 'continuous') xx[j] <- datj
          else {
            levj <- levels(datj)
            if(! length(levj)) levj <- unique(datj)
            Lev <- c(Lev, levj)
            xx[j] <- as.character(datj)
          }
        }
        if(type == 'discrete') {
          Lev <- unique(Lev)
          xx  <- factor(xx, Lev)
        }
        dat$.xx. <- xx
        if(length(groups)) dat$.co. <- as.factor(dat[[groups]])

        ylimc <- limc(ylim., 'y')
        if(type == 'continuous') {
          if(length(groups)) g <- 
            sprintf('ggplot(dat, aes(x=.xx., y=yhat, %s=%s)) +
                     labs(x=NULL, y=%s) + %s',
                    aestype[1], groups[1], expch(ylab), ylimc)
          else
            g <- sprintf("ggplot(dat, aes(x=.xx., y=yhat)) +
                         labs(x=NULL, y=%s) + %s",
                         expch(ylab), ylimc)
          
          g <- c(g, if(length(layout))
                      sprintf("facet_wrap(~ .predictor., scales='free_x',
                                ncol=%s)",
                              layout[2]) else
                 "facet_wrap(~ .predictor., scales='free_x')", "geom_line()")
          if(conf.int) {
            h <- 
              if(conf == 'fill')
                sprintf("geom_ribbon(aes(x=.xx., ymin=lower, ymax=upper),%s)",
                        ribbonargs)
              else c("geom_line(aes(x=.xx., y=lower))",
                     "geom_line(aes(x=.xx., y=upper))")
            g <- c(g, h)
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
            g <- c(g,
                   sprintf("dohist(%s, predictions=dat, data=rdata, ylim=%s)",
                           form, deparse(ylim.)))
          }
        } else {   # discrete x
          if(length(groups)) g <- 
            c(sprintf('ggplot(dat, aes(x=yhat, y=.xx., %s=%s))',
                      aestype[1], groups[1]),
              sprintf("labs(x=%s, y=NULL)", expch(ylab)))
          else
            g <- c("ggplot(dat, aes(x=yhat, y=.xx.))",
                   sprintf("labs(x=%s, y=NULL)", expch(ylab)))
          if(! maddlayer) g <- c(g, addlayer)
          g <- c(g, limc(ylim., 'x'),
                 "facet_wrap(~ .predictor., scales='free_y')",
                 "geom_point()")
          if(conf.int) g <- c(g,
              "geom_errorbarh(aes(y=.xx., xmin=lower, xmax=upper), height=0)")
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
            datj <- data[j,, drop=FALSE]
            xv <- datj[, iv]
            xx <- switch(type, continuous = xv,
                               discrete   = as.numeric(xv))
            yy <- datj[, 'yhat']
            if(conf.int) {
              xx <- c(xx, xx, xx)
              yy <- c(yy, datj[, 'lower'], datj[, 'upper'])
            }
            xlim. <- if(is.factor(xv)) c(1, length(levels(xv)))
             else range(pretty(xv))
            tan <- tanova(iv, xx, yy, xlim., ylim., dataOnly=TRUE,
                          flip=type=='discrete', empty=type == 'discrete')
            .xx. <- c(.xx., tan$x)
            yhat <- c(yhat, tan$y)
            .label. <- c(.label., tan$label)
            hjust <- c(hjust, tan$hjust)
            vjust <- c(vjust, tan$vjust)
          }
          .anova. <- data.frame(.predictor.=v, .xx., yhat, .label.,
                                hjust, vjust)
          g <- c(g, "geom_text(aes(label=.label., hjust=hjust, vjust=vjust),
                             size=size.anova,
                             data=.anova., parse=TRUE, show_guide=FALSE)")
        }
        g <- paste(g, collapse=' + ')
        if(ggexpr) return(g)
        g <- eval(parse(text = g))
        if(vnames == 'labels') g <- facet_wrap_labeller(g, pmlabel[v])
        g
      }
      
      gcont <- if(any(! isdis)) dogroup('continuous')
      gdis  <- if(any(  isdis)) dogroup('discrete')
      if(ggexpr) return(list(continuous=gcont, discrete=gdis))
      r <- mean(! isdis)
      return(if(length(gcont) && length(gdis))
               switch(sepdiscrete,
                      list = list(continuous=gcont, discrete=gdis),
                      vertical = arrGrob(gcont, gdis,         heights=c(r, 1-r)),
                      horizontal=arrGrob(gcont, gdis, nrow=1, widths =c(r, 1-r)))
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
      xlim. <- if(ll) c(1, ll) else range(pretty(z))
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
        g <- sprintf(
          'ggplot(zz, aes(x=.xx., y=.yhat, %s=.cond))', aestype[1])
      }
      else g <- 'ggplot(zz, aes(x=.xx., y=.yhat))'
      
      xdiscrete <- is.factor(z) || is.character(z) ||
        length(unique(z[!is.na(z)])) <= nlevels
      flipped <- FALSE
      if(xdiscrete) {
        if(flipxdiscrete && ! is.numeric(z)) {
          g <- c(g, 'coord_flip()')
          flipped <- TRUE
        }
        g <- c(g, 'geom_point()',
               if(length(type) && type %in% c('l', 'b')) 'geom_line()')
        if(is.numeric(z))
          g <- c(g, sprintf('scale_x_discrete(breaks=%s)', deparse(unique(z))))
      } else {
        if(length(type))
          g <- c(g, switch(type, p='geom_point()', l='geom_line()',
                                 b='geom_point() + geom_line()'))
        else g <- c(g, 'geom_line()')
      }

      ## Need the following or geom_ribbon will improperly clip regions
      if(flipped) g <- c(g, limc(ylim., 'y')) else
      g <- c(g, sprintf('coord_cartesian(ylim=%s)', deparse(ylim.)))
      g <- c(g, sprintf('labs(x=%s, y=%s)',
                        expch(xl), expch(ylab)),
             "theme(plot.margin = grid::unit(rep(0, 4), 'cm'))")
      ## use rep(.1, 4) if using print(..., viewport=...) for multiple plots
      if(length(groups)) {
        # if(nr == 1 && nc == 1) {
        if(jplot == 1) {
          colFun <- if(aestype[1] == 'color') colorscale else
           get(paste('scale', aestype[1], 'discrete', sep='_'))
          groupLabel <- glabel(groups[1])
          g <- c(g, if(aestype[1] == 'size')
                      sprintf("colFun(name=%s, range=c(.2, 1.5))",
                              expch(groupLabel)) else
                 sprintf("colFun(name=%s)", expch(groupLabel)))
          g <- c(g, sprintf("theme(legend.position='%s')",
                            legend.position))
          
        } else g <- c(g, "theme(legend.position='none')")
      }
      xa <- if(conf.int) c(zz$.xx., zz$.xx., zz$.xx.) else zz$.xx.
      ya <- if(conf.int) c(zz$.yhat, zz$lower, zz$upper) else zz$.yhat
      g  <- c(g,
              sprintf("tanova(w, xa, ya, %s, %s, flip=FALSE)",
                      deparse(xlim.), deparse(ylim.)))   ## was flip=flipped
      
      if(! maddlayer) g <- c(g, addlayer)
      if(conf.int) {
        h <- 
          if(ll || xdiscrete)
            "geom_errorbar(aes(ymin=lower, ymax=upper), width=0)"
          else {
            if(conf == 'fill')
              sprintf("geom_ribbon(aes(ymin=lower, ymax=upper), %s)",
                      ribbonargs)
            else c("geom_line(aes(x=.xx., y=lower))",
                   "geom_line(aes(x=.xx., y=upper))")
            ## geom_ribbon with fill=NA draws vertical lines at
            ## ends of confidence regions
          }
        g <- c(g, h)
      }
      
      if(! missing(formula))
        g <- c(g, sprintf("facet_grid(%s)", deparse(formula)))
      
      if(! is.factor(z) && length(rdata) && w %in% names(rdata)) {
        rdata$.xx. <- rdata[[w]]
        if(length(.co)) {
          rdata$.cond <- rdata[[groups]]
          form <- '.yhat ~ .xx. + .cond'
        } else form <- '.yhat ~ .xx.'
        g <- c(g, sprintf("dohist(%s, predictions=zz, data=rdata, ylim=%s)",
                          form, deparse(ylim.)))
      }
      # print(g, vp = viewport(layout.pos.row=nr, layout.pos.col=nc))
      g <- paste(g, collapse = ' + ')
      if(ggexpr) return(g)
      g <- eval(parse(text=g))
      Plt[[jplot]] <- g
    }
    Plt <- if(jplot == 1) Plt[[1]]
     else
       do.call(arrGrob, c(Plt, list(ncol=layout[2])))
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
      j <- if(! length(groups)) perim(xv, NULL) else
                                perim(xv, data[[groups[1]]])
      data$yhat[! j] <- NA
      if(conf.int) data$lower[! j] <- data$upper[! j] <- NA
    }
    ae <- paste('aes(x=', xn, ', y=yhat', sep='')
    if(length(groups)) for(j in 1 : length(groups))
      ae <- paste(ae, ', ', aestype[j], '=', groups[j], sep='')
    ae <- eval(parse(text=paste(ae, ')', sep='')))
    g <- c("ggplot(data, ae)", sprintf("labs(x=%s, y=%s)",
                                       expch(xlab), expch(ylab)))

    flipped <- FALSE
    if(xdiscrete) {
      if(flipxdiscrete && ! is.numeric(xv)) {
        g <- c(g, "coord_flip()")
        flipped <- TRUE
      }
      g <- c(g, "geom_point()", 
             if(length(type) && type %in% c('l', 'b')) "geom_line()" )

      if(conf.int) g <- c(g,
                          "geom_errorbar(aes(ymin=lower, ymax=upper), width=0)")
    } else {
      if(length(type))
        g <- c(g, switch(type, p="geom_point()", l="geom_line()",
                               b="geom_point() + geom_line()"))
      else g <- c(g, "geom_line()")

      if(length(groups)) {
        for(j in 1 : length(groups)) {
          colFun <- if(aestype[j] == 'color') colorscale else
           get(paste('scale', aestype[j], 'discrete', sep='_'))
          groupLabel <- glabel(groups[j], j)
          g <- c(g, if(aestype[j] == 'size')
                 sprintf("colFun(name=%s, range=c(.2, 1.5))",
                         expch(groupLabel)) else
                 sprintf("colFun(name=%s)", expch(groupLabel)))
        }
        g <- c(g, sprintf("theme(legend.position='%s')",
                          legend.position))
      }
      if(conf.int) {
        h <- 
          if(conf == 'fill')
            sprintf("geom_ribbon(data=data, aes(ymin=lower, ymax=upper),%s)",
                    ribbonargs)
          else c(sprintf('geom_line(data=data, aes(x=%s, y=lower))', xn),
                 sprintf('geom_line(data=data, aes(x=%s, y=upper))', xn))
        g <- c(g, h)
      }  # end if(conf.int)
    }
    if(! maddlayer) g <- c(g, addlayer)

    g <- c(g, if(flipped) sprintf("ylim(%s)", deparse(ylim.)) else
           sprintf("coord_cartesian(ylim=%s)", deparse(ylim.)))

    xa <- if(conf.int) c(xv, xv, xv) else xv
    ya <- if(conf.int) c(data$yhat, data$lower, data$upper) else data$yhat
    if(missing(xlim.))
      xlim. <- if(is.factor(xv)) c(1 , length(levels(xv)))
       else
         range(pretty(xv))

    g <- c(g, sprintf("tanova(xn, xa, ya, %s, %s, flip=FALSE)",
                      deparse(xlim.), deparse(ylim.)))  # was flip=flipped
                    
    if(! is.factor(xv) && length(rdata) && xn %in% names(rdata)) {
      form <- paste('yhat', xn, sep='~')
      if(length(groups)) form <- paste(form, groups[1], sep='+')
      g <- c(g, sprintf("dohist(%s, predictions=data, data=rdata, ylim=%s)",
                        form, deparse(ylim.)))
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
      } else formula <- deparse(formula)
      g <- c(g, sprintf("facet_grid(%s)", formula))
    }
    g <- paste(g, collapse=' + ')
    if(ggexpr) return(g)
    g <- eval(parse(text=g))
    if(length(sub)) g <- footnote(g, sub)
    g
  }
}

utils::globalVariables(c('.xx.', '.yhat', 'lower', 'upper', 'groupLabel'))
