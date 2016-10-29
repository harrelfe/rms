ggplot.Predict <-
  function(data, mapping, formula=NULL, groups=NULL,
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
           type=NULL, ggexpr=FALSE, height=NULL, width=NULL, ..., environment)
{
  isbase <- Hmisc::grType() == 'base'   ## vs. 'plotly'
  if(! isbase && length(anova))
    stop('anova not yet implemented for grType plotly')
  lhw <- length(height) + length(width)
  if(isbase && lhw)
    warning('height and width ignored for non-plotly graphics')

  plrend <- if(isbase) function(obj, ...) obj
            else
              function(obj, final=TRUE) {
                if(final && (length(width) > 0 || length(height) > 0))
                  plotly::ggplotly(obj, height=height, width=width)
                else
                  plotly::ggplotly(obj)
              }

  comb <- function(plist, nrow=1, ncol=1, ...) {
    ## Note: subplot does not take an ncols argument
    if(isbase) return(do.call(arrGrob,
                              c(plist, list(nrow=nrow, ncol=ncol), list(...))))
    z <- 
      do.call(plotly::subplot,
              c(plist, list(nrows=nrow,
                            titleX=TRUE, titleY=TRUE, margin=0.065),
                list(...)))
    if(lhw) {
      z <- plotly::plotly_build(z)
      z$layout$height <- height
      z$layout$width  <- width
    }
    # if(lhw) z <- layout(z, height=height, width=width)  # also works
    z
  }
  
  if(! length(formula) && ! missing(mapping)) formula <- mapping
  ## .xlim, .ylim instead of xlim, ylim to distinguish from ggplot functions
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

  ribbonargs <- sprintf("alpha=0.2, linetype=0, fill=I('%s'),show.legend=FALSE",
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
    pmlabel[i] <-
      if(isbase) as.character(labelPlotmath(label[i], units[i]))
      else markupSpecs$html$varlabel(label[i], units[i])
  
  if(predpres)
    data$.Predictor. <- if(vnames != 'labels') data$.predictor.
                        else pmlabel[as.character(data$.predictor.)]
    
  glabel <- function(gname, j=1, chr=FALSE) {
    r <-  
      if(! length(legend.label))
        if(isbase) parse(text=pmlabel[gname]) else pmlabel[gname]
      else if(is.logical(legend.label)) ''
      else legend.label[j]
    if(is.expression(r)) {
      if(chr) r <- sprintf('expression(%s)', as.character(r))
    } else {
      qc <- if(length(grep("'", r))) '"' else "'"
      r <- paste0(qc, r, qc)
      }
    r
  }

  ## Function to create expression( ) or "" depending on argument
  expch <- function(x, chr=FALSE) {
    if(! length(x)) 'NULL' else
      if(is.expression(x)) {
        if(chr) sprintf('expression(%s)', as.character(x)) else x
      } else if(grepl('expression\\(', x)) x else deparse(x)
    }

  ## Function to construct xlim() or ylim() call
  limc <- function(limits, which)
    sprintf("%slim(%s, %s)", which, limits[1], limits[2])

  xlimc <- if(missing(xlim.)) ''
           else
             paste('+', limc(xlim., 'x'))
  
  if(! missing(subset)) {
    subset <- eval(substitute(subset), data)
    data   <- data[subset,, drop=FALSE]
  }

  if(length(groups) == 1 && is.logical(groups) && ! groups) groups <- NULL
  else if(length(groups)) {
    if(length(groups) > 2 || !is.character(groups) ||
       any(groups %nin% names(data)))
      stop('groups must be one or two predictor names')
    ## geom_ribbon will not handle two aesthetics
    if(length(groups) == 2) conf.int <- FALSE
  } else if(! predpres && length(varying) > 1) groups <- varying[2]

  ## Make all grouping variables discrete for proper aesthetic mapping
  if(length(groups)) for(v in groups) data[[v]] <- as.factor(data[[v]])
  
  if(missing(ylab))
    ylab <- if(isbase) info$ylabPlotmath else info$ylabhtml
  if(! length(data$lower)) conf.int <- FALSE
  
  if(missing(ylim.))
    ylim. <- range(pretty(
      if(conf.int) c(data$yhat, data$lower, data$upper)
      else data$yhat), na.rm=TRUE)

  if(missing(adj.subtitle)) adj.subtitle <- length(adjust) > 0
  sub <- if(adj.subtitle && length(adjust)==1)
           paste0('Adjusted to:', adjust) else NULL
  cap <- expch(sub, chr=TRUE)
  ## Won't need the following when ggplot2 really gets subtitle captions working
  thm <- if(length(sub)) 'theme(plot.title=element_text(size=8, hjust=1))'
  ## hjust seems to be ignored some of the time
  
  ## ggplot2 is supposed to implement labs(subtitle= caption=) but
  ## neither of these work as of 2016-07-18.  labs(title=) used for now.

  tanova <- if(length(anova))
    function(name, x, y, xlim, ylim, flip=FALSE,
             empty=FALSE, dataOnly=FALSE)
      annotateAnova(name, plotmathAnova(anova, pval),
                    x, y, ggplot=TRUE, xlim=xlim, ylim=ylim,
                    size=size.anova, flip=flip, empty=empty, dataOnly=dataOnly)
    else function(...) {}

  ## See http://bigdata-analyst.com/best-way-to-add-a-footnote-to-a-plot-created-with-ggplot2.html
  ## size is in mm
#  footnote <- function(object, text, size=2.5, color=grey(.5))
#    arrGrob(object, sub = grid::textGrob(text, x = 1, hjust = 1.01,
#          vjust=0.1, gp = grid::gpar(fontsize =size/0.3527778 )))
  
  if(predpres) {   ## User did not specify which predictors to plot; all plotted
    data$.predictor.  <- factor(data$.predictor.)

    if(sepdiscrete != 'no') {
      ## From http://stackoverflow.com/questions/11979017
      ## Changed to assume that each element of labels is a character string
      ## of the form "expression(....)"
      if(FALSE) facet_wrap_labeller <- function(gg.plot, labels=NULL) {
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
        xx <- switch(type,
                     continuous = numeric(  nrow(dat)),
                     discrete   = character(nrow(dat)) )

        lbr <- if(! isbase || vnames != 'labels') ''
               else
                 ', labeller=label_parsed'
        
        ## Prepare to create a "super factor" variable by concatenating
        ## all levels of all categorical variables keeping original orders
        ## firstLev is first level for each discrete predictor
        ## Thought was needed with anova but geom_text will take a numeric
        ## x or y coordinate where factor levels seem to be spaced at 1.0
        Lev <- character()
        ## firstLev <- character(length(v))
        ## names(firstLev) <- v

        for(iv in v) {
          j <- which(p == iv)
          datj <- dat[j, iv]
          if(type == 'continuous') {
            xx[j] <- datj
            ## firstLev[iv] <- ''
          }
          else {
            levj <- levels(datj)
            if(! length(levj)) levj <- unique(datj)
            Lev <- c(Lev, levj)
            xx[j] <- as.character(datj)
            ## firstLev[iv] <- levj[1]
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
          if(length(groups))
            g <- 
              sprintf('ggplot(dat, aes(x=.xx., y=yhat, %s=%s)) +
                     labs(x=NULL, y=%s, title=%s) + %s %s',
                     aestype[1], groups[1], expch(ylab, chr=TRUE), cap,
                     ylimc, xlimc)
          else
            g <- sprintf("ggplot(dat, aes(x=.xx., y=yhat)) +
                         labs(x=NULL, y=%s, title=%s) + %s %s",
                         expch(ylab, chr=TRUE), cap, ylimc, xlimc)
          
          g <- c(g, if(length(layout))
                      sprintf("facet_wrap(~ .Predictor., scales='free_x',
                                ncol=%s%s)",
                              layout[2], lbr)
                    else
                      sprintf("facet_wrap(~ .Predictor., scales='free_x'%s)",
                              lbr), "geom_line()")
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
                             timevar='.Predictor.',
                             varying=rv, times=rv)
            rdata$.Predictor. <- pmlabel[rdata$.Predictor.]
            form <- 'yhat ~ .xx. + .Predictor.'
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
              sprintf("labs(x=%s, y=NULL, title=%s)",
                      expch(ylab, chr=TRUE), cap))
          else
            g <- c("ggplot(dat, aes(x=yhat, y=.xx.))",
                   sprintf("labs(x=%s, y=NULL, title=%s)",
                           expch(ylab, chr=TRUE), cap))
          if(! maddlayer) g <- c(g, addlayer)
          g <- c(g, limc(ylim., 'x'),
                 sprintf("facet_wrap(~ .Predictor., scales='free_y'%s)", lbr),
                 "geom_point()")
          if(conf.int) g <- c(g,
              "geom_errorbarh(aes(y=.xx., xmin=lower, xmax=upper), height=0)")
        }
        ## anova annotations need to be created for all variables being
        ## plotted with faceting, and annotation information must be
        ## based on a dataset with the information and the .Predictor.
        ## variable, and geom_text() must be used instead of annotate()
        ## See http://stackoverflow.com/questions/2417623
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
            ## .xx. <- c(.xx., if(type == 'discrete') firstLev[iv] else tan$x)
            .xx. <- c(.xx., tan$x)
            yhat <- c(yhat, tan$y)
            .label. <- c(.label., tan$label)
            hjust <- c(hjust, tan$hjust)
            vjust <- c(vjust, tan$vjust)
          }
          .anova. <-
            data.frame(.Predictor. = if(vnames != 'labels') v else pmlabel[v],
                       .xx., yhat, .label.,
                       hjust, vjust)
          g <- c(g, sprintf("geom_text(aes(label=.label., hjust=hjust, vjust=vjust),
                             size=size.anova, nudge_y=%s,
                             data=.anova., parse=TRUE, show.legend=FALSE)",
                            if(type == 'discrete') -0.25 else 0))
        }
        g <- c(g, thm)
        g <- paste(g, collapse=' + ')
        if(ggexpr) return(g)
        g <- eval(parse(text = g))
        g
      }    # end dogroup function
      
      gcont <- if(any(! isdis)) dogroup('continuous')
      gdis  <- if(any(  isdis)) dogroup('discrete')
      if(ggexpr) return(list(continuous=gcont, discrete=gdis))
      r <- mean(! isdis)

      return(if(length(gcont) && length(gdis))
               switch(sepdiscrete,
                      list = list(continuous = plrend(gcont),
                                  discrete   = plrend(gdis )),
                      vertical  = comb(list(plrend(gcont, final=FALSE),
                                            plrend(gdis,  final=FALSE)),
                                       nrow=2, heights=c(r, 1-r)),
             horizontal= comb(list(plrend(gcont, final=FALSE),
                                   plrend(gdis,  final=FALSE)),
                              ncol=2, widths =c(r, 1-r)))
             else if(length(gcont)) plrend(gcont) else plrend(gdis))
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
        else if(np <=20) c(4,5)
        else ceil(rep(sqrt(np), 2))        
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
      xl <- if(vnames == 'labels') {
              if(isbase) parse(text=pmlabel[w])
              else pmlabel[w]
            }
            else w
      zz <- data.frame(.xx.=z, .yhat=yhat)
      if(length(formula))
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
      g <- c(g, sprintf('labs(x=%s, y=%s, title=%s)',
                        expch(xl, chr=TRUE), expch(ylab, chr=TRUE), cap),
             "theme(plot.margin = grid::unit(rep(0, 4), 'cm'))")
      ## use rep(.1, 4) if using print(..., viewport=...) for multiple plots
      if(length(groups)) {
#### ??
        # if(nr == 1 && nc == 1) {
        if(jplot == 1) {
          colFun <- if(aestype[1] == 'color') colorscale else
           get(paste('scale', aestype[1], 'discrete', sep='_'))
          groupLabel <- glabel(groups[1], chr=TRUE)
          g <- c(g, if(aestype[1] == 'size')
                      sprintf("colFun(name=%s, range=c(.2, 1.5))",
                              groupLabel) else
                 sprintf("colFun(name=%s)", groupLabel))
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
      
      if(length(formula))
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
      g <- c(g, thm)
      g <- paste(g, collapse = ' + ')
      if(ggexpr) return(g)
      g <- eval(parse(text=g))
      Plt[[jplot]] <- g
    }

res <- if(jplot == 1) plrend(Plt[[1]])
       else {
             for(j in 1 : jplot) Plt[[j]] <- plrend(Plt[[j]], final=FALSE)
             comb(Plt, nrow=layout[1], ncol=layout[2])
           }
#    if(length(sub)) {
#      Plt <- if(isbase) footnote(Plt, sub, size=size.adj)
#             else
#               plotly::layout(p, title=sub, margin=0.03)
#      }
    return(res)
    
  } else  { # .predictor. not included; user specified predictors to show
    v  <- varying
    xn <- v[1]    ## name of x-axis variable (first variable given to Predict)
    if(missing(xlab))
      xlab <- if(isbase) parse(text = pmlabel[xn]) else pmlabel[xn]
    xv <- data[[xn]]
    xdiscrete <- is.factor(xv) || is.character(xv) ||
      length(unique(xv[!is.na(xv)])) <= nlevels
    if(length(perim)) {
      j <- if(! length(groups)) perim(xv, NULL) else
                                perim(xv, data[[groups[1]]])
      data$yhat[! j] <- NA
      if(conf.int) data$lower[! j] <- data$upper[! j] <- NA
    }
    ae <- paste0('aes(x=', xn, ', y=yhat')
    if(length(groups)) for(j in 1 : length(groups))
      ae <- paste0(ae, ', ', aestype[j], '=', groups[j])
####    ae <- eval(parse(text=paste0(ae, ')')))
    ae <- paste0(ae, ')')
    g <- c(sprintf("ggplot(data, %s)", ae),
           sprintf("labs(x=%s, y=%s, title=%s) %s",
              expch(xlab, chr=TRUE), expch(ylab, chr=TRUE), cap, xlimc))

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
#          colFun <- if(aestype[j] == 'color') colorscale else
#           get(paste('scale', aestype[j], 'discrete', sep='_'))
          colFun <- if(aestype[j] == 'color') 'colorscale'
                    else
                      paste('scale', aestype[j], 'discrete', sep='_')
          groupLabel <- glabel(groups[j], j, chr=TRUE)
#??          g <- c(g, if(aestype[j] == 'size')
#                 sprintf("colFun(name=%s, range=c(.2, 1.5))",
#                         groupLabel) else
#                 sprintf("colFun(name=%s)", groupLabel))
          g <- c(g, if(aestype[j] == 'size')
                      sprintf('%s(name=%s, range=c(.2, 1.5))',
                              colFun, groupLabel)
                    else
                      sprintf('%s(name=%s)', colFun, groupLabel))
          
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
      if(! length(formula)) {
        k <- length(f)
        formula <- if(k == 1) paste('~', f[1])
                   else if(k == 2) paste(f[1], f[2], sep='~')
                   else if(k == 3) paste(f[1], '~', f[2], '*', f[3])
                   else if(k == 4) paste(f[1], '*', f[2], '~', f[3], '*', f[4])
                   else stop('too many varying variables to use faceting')
      } else formula <- deparse(formula)
      g <- c(g, sprintf("facet_grid(%s)", formula))
    }
    g <- c(g, thm)
    g <- paste(g, collapse=' + ')
    if(ggexpr) return(g)
    g <- plrend(eval(parse(text=g)))
#    if(length(sub)) g <- if(isbase) footnote(g, sub)
#                         else
#                           plotly::layout(g, title=sub, margin=0.03)
    ## Could not get layout(g, annotations=...) to work
    return(g)
  }
}

utils::globalVariables(c('.xx.', '.yhat', 'lower', 'upper', 'groupLabel'))
