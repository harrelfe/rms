survplotp <- function(fit, ...) UseMethod("survplotp")

survplotp.npsurv <-
  function(fit, xlim, 
           ylim, xlab, ylab, time.inc, state=NULL,
           conf=c("bands", "none"), mylim=NULL, abbrev.label=FALSE,
           col=colorspace::rainbow_hcl,
           levels.only=TRUE,
           loglog=FALSE, fun=function(y) y, aehaz=FALSE, times=NULL,
           logt=FALSE, pr=FALSE, ...) {
    
    conf     <- match.arg(conf)
    conf.int <- fit$conf.int
    if(!length(conf.int) | conf == "none") conf.int <- 0

    if(loglog)
      fun <- function(y) logb(-logb(ifelse(y == 0 | y == 1, NA, y)))

    mstate <- inherits(fit, 'survfitms')
    if(mstate) fun <- function(y) 1 - y
    ## Multi-state model for competing risks

    z <- seq(.1, .9, by = .1)
    funtype <- if(loglog) 'loglog'
               else if(all(fun(z) == z)) 'identity'
               else if(all(abs(fun(z) - (1 - z)) < 1e-6)) 'inverse'
               else if(loglog) 'loglog'
               else 'other'

    cylim <- function(ylim)
      if(length(mylim)) c(min(ylim[1], mylim[1]), max(ylim[2], mylim[2]))
      else ylim

    mu      <- markupSpecs$html
    nbsp    <- htmlSpecial('nbsp')

    if(funtype %in% c('identity', 'inverse'))
      survdiffplotp <-
        function(fit, fun, xlim, 
                 conf.int, convert=function(f) f, pobj)
      {
        if(length(fit$strata) != 2)
          stop('must have exactly two strata')
      
      h <- function(level, f) {
        i <- f$strata == levels(f$strata)[level]
        tim   <- f$time[i]
        surv  <- fun(f$surv[i])
        se    <- f$std.err[i]
        list(time=tim, surv=surv, se=se)
      }
      
      times <- sort(c(0, unique(fit$time)))
      times <- times[times >= xlim[1] & times <= xlim[2]]
      
      f <- convert(summary(fit, times=times, print.it=FALSE, extend=TRUE))
      a <- h(1, f)
      b <- h(2, f)
      
      if(! identical(a$time, b$time)) stop('program logic error')
      time  <- a$time
      surv  <- (a$surv + b$surv) / 2
      se    <- sqrt(a$se^2 + b$se^2)
      
      z  <- qnorm((1 + conf.int) / 2)
      lo <- pmax(0, surv - 0.5 * z * se)
      hi <- pmin(1, surv + 0.5 * z * se)
      k <- ! is.na(time + lo + hi)
      list(times=time[k], lower=lo[k], upper=hi[k])
    }
    
    
    fit.orig <- fit
    units <- fit$units
    if(!length(units)) units <- "Day"
    maxtime <- fit$maxtime
    if(! length(maxtime)) maxtime <- max(fit$time)
    mintime <- min(fit$time, 0)
    pret    <- pretty(c(mintime, maxtime))
    maxtime <- max(pret)
    mintime <- min(pret)
    if(missing(time.inc)) {
      time.inc <- switch(units, Day=30, Month=1, Year=1,
                                (maxtime - mintime) / 10)
      if(time.inc > maxtime) time.inc <- (maxtime - mintime) / 10
    }

    if(mstate) {
      ## Multi-state model for competing risks
      if(missing(state))
        stop('state must be given when response is a multi-state/competing risk object from Surv()')
      if(length(state) != 1) stop('at present state can only be a single state')
      states <- fit$states
      if(state %nin% states) stop(paste('state is not in',
                                        paste(states, collapse=', ')))
    }

    if(missing(ylab))
      ylab <- switch(funtype,
                     identity = 'Survival Probability',
                     inverse  = if(mstate) paste('Cumulative Incidence of',
                                                 upFirst(state))
                                else 'Cumulative Incidence',
                     loglog   = 'log(-log Survival Probability)',
                     other    = '')

    un <- fit$units
    if(un != '') un <- paste0(un, 's')
    if(missing(xlab))
      xlab <- if(logt) paste0("log Follow-up Time in ", un)
              else 
                mu$varlabel('Follow-up Time', un)
  
    if(missing(xlim)) 
      xlim <- if(logt) logb(c(maxtime / 100, maxtime)) else c(mintime, maxtime)

    convert <- if(mstate) {
      istate    <- match(state, states)
      conv <- function(f, istate) {
        f$surv    <- 1 - f$pstate [, istate]
        f$lower   <- 1 - f$lower  [, istate]
        f$upper   <- 1 - f$upper  [, istate]
        f$std.err <-     f$std.err[, istate]
        f$n.risk  <-     f$n.risk[, ncol(f$n.risk)]
        if(all(f$n.risk == 0))
          stop('expected n.risk to be last column of n.risk matrix for competing risks')
        f
      }
      formals(conv) <- list(f=NULL, istate=istate)
      conv
               }
               else function(f) f
    
    fit <- convert(fit)
  
    origsurv <- fit$surv
    if(funtype != 'identity') {
      fit$surv <- fun(fit$surv)
      fit$surv[is.infinite(fit$surv)] <- NA
      ##  handle e.g. logit(1) - Inf would mess up ylim in plot()
      if(conf.int > 0) {
        fit$lower <- fun(fit$lower)
        fit$upper <- fun(fit$upper)
        fit$lower[is.infinite(fit$lower)] <- NA
        fit$upper[is.infinite(fit$upper)] <- NA
        if(missing(ylim))
          ylim <- cylim(range(c(fit$lower, fit$upper), na.rm=TRUE))
      }
      else if(missing(ylim)) ylim <- cylim(range(fit$surv, na.rm=TRUE))
    }
    else if(missing(ylim)) ylim <- c(0, 1)

    olev <- slev <- names(fit$strata)
    if(levels.only) slev <- gsub('.*=', '', slev)
    sleva <- if(abbrev.label) abbreviate(slev) else slev
    ns <- length(slev)
    slevp <- ns > 0
    
    ns  <- max(ns, 1)
    
    if(is.function(col)) col <- col(ns)
    
    y <- 1 : ns
    strat <- if(ns == 1) rep(1, length(fit$time))
             else rep(1 : ns, fit$strata)
    
    stime <- sort(unique(c(0, fit.orig$time)))
    stime <- stime[stime >= mintime & stime <= maxtime]
#    v <- convert(summary(fit.orig, times=stime, print.it=FALSE))
#    vs <- if(ns > 1) as.character(v$strata)
    ## survival:::summary.survfit was not preserving order of strata levels

    nevents <- totaltime <- numeric(ns)
    cuminc  <- character(ns)

    p <- plotly::plot_ly(...)

    pl <- function(p, x, y, n.risk=NULL, col, slev, type='est') {
      sname  <- if(ns == 1) '' else slev
      snames <- if(sname == '') '' else paste0(sname, ' ')
      d <- paste0('Difference<br>', mu$half(), ' ',
                  conf.int, ' CL')
      nam   <- switch(type,
                      est   = sname,
                      lower = paste0(snames, conf.int, ' CL'),
                      upper = paste0(snames, conf.int, ' CL', type),
                      'diff lower' = d,
                      'diff upper' = paste0(d, type))
      lg <- switch(type,
                   est = 'Estimates',
                   lower = paste0(snames, 'CL'),
                   upper = paste0(snames, 'CL'),
                   'diff lower' = 'Difference',
                   'diff upper' = 'Difference')
      rx <- format(round(x, 3))
      ry <- format(round(y, 3))
      txt <- switch(type,
                    est   = paste0('t=', rx,
                                   '<br>Probability=', ry,
                                   if(length(n.risk)) '<br>At risk:', n.risk),
                    lower = paste0('t=', rx, '<br>Lower:', ry),
                    upper = paste0('t=', rx, '<br>Upper:', ry),
                    'diff lower' = NULL,
                    'diff upper' = NULL)

      ## Mark in text the point estimates that correspond to every time.inc
      if(type == 'est' && ! logt) {
        nicet <- seq(mintime, maxtime, by=time.inc)
        nicet <- nicet[nicet > 0]
        for(ti in nicet) {
          if(any(abs(ti - x) < 1e-6)) next  # nice time already covered
          k <- which(x < ti); k <- max(k)
          txt[k] <- paste0(txt[k], '<br>(Also for t=', ti, ')')
        }
      }
      
      fcol <- plotly::toRGB(col, 0.2)
      vis  <- if(ns == 2 && type %in% c('lower', 'upper'))
                'legendonly' else TRUE
      ln <- if(type == 'est') list(shape='hv', color=col)
            else list(shape='hv', color=col, width=0)

      dat <- if(length(txt)) data.frame(x, y, txt) else data.frame(x, y)
      up <- type %in% c('upper', 'diff upper')
      p <- if(length(txt))
             plotly::add_lines(p, x=~ x, y=~ y, text=~ txt, data=dat,
                               hoverinfo='text', line=ln,
                               fillcolor=fcol,
                               fill=if(type %in% c('upper', 'diff upper'))
                                      'tonexty' else 'none',
                               visible=vis, legendgroup=lg,
                               name=nam, showlegend=! up)
           else
             plotly::add_lines(p, x=~ x, y=~ y, data=dat,
                               hoverinfo='none', line=ln,
                               fillcolor=fcol,
                               fill=if(type %in% c('upper', 'diff upper'))
                                      'tonexty' else 'none',
                               visible=vis, legendgroup=lg,
                               name=nam, showlegend=! up)

    }
                        
    for(i in 1 : ns) {
      st <- strat == i
      time         <- fit$time[st]
      surv         <- fit$surv[st]
      lower        <- fit$lower[st]
      upper        <- fit$upper[st]
      osurv        <- origsurv[st]
      n.risk       <- fit$n.risk[st]

      if(! logt && xlim[1] ==0 && all(time > xlim[1])) {
        time   <- c(xlim[1], time)
        surv   <- c(fun(1), surv)
        lower  <- c(fun(1), lower)
        upper  <- c(fun(1), upper)
        osurv  <- c(1, osurv)
        n.risk <- c(fit$n[i], n.risk)
      }
      
      ## nevents[i]   <- sum(fit$n.event[st])
      ## nrsk         <- fit$n.risk[st]
      ## neachtime    <- c(- diff(nrsk), min(nrsk))
      ## totaltime[i] <- sum(neachtime * time)

      nevents[i] <-
        if(mstate) {
          if(ns == 1) fit$numevents[, state]
          else fit$numevents[olev[i], state]
        } else {
          if(ns == 1) fit$numevents else fit$numevents[olev[i]]
        }
      totaltime[i] <- if(ns == 1) fit$exposure else fit$exposure[olev[i]]
      if(length(times)) {
        cumi <- 1. - approx(time, osurv, xout=times, method='constant')$y
        noun <- units %in% c('', ' ')
        cuminc[i]   <- paste(round(cumi, 3), collapse=', ')
      }
      if(logt) time <- logb(time)
      ##don't let step function go beyond x-axis -
      ##this cuts it off but allows step to proceed axis end
      if(max(time) > xlim[2]) {
        srvl <- surv[time <= xlim[2] + 1e-6]
        s.last <- srvl[length(srvl)]
        k <- time < xlim[2]
        time <- c(time[k], xlim[2])
        surv <- c(surv[k], s.last)
        n.risk <- c(n.risk[k], n.risk[length(srvl)])
        if(conf.int > 0) {
          low.last <- lower[time <= xlim[2] + 1e-6]
          low.last <- low.last[length(low.last)]
          up.last  <- upper[time <= xlim[2] + 1e-6]
          up.last  <- up.last[length(up.last)]
          lower   <- c(lower[k], low.last)
          upper   <- c(upper[k], up.last)
        }
      }

      if(logt) p <- pl(p, time, surv, n.risk, col=col[i], slev=sleva[i])
      else {
        xxx <- time
        yyy <- surv
        nr  <- n.risk
        if(mintime < min(time)) {
          xxx <- c(mintime, time)
          yyy <- c(fun(1), surv)
          nr  <- c(fit$n[i], n.risk)
        }
        p <- pl(p, xxx, yyy, nr, col=col[i], slev=sleva[i])
      }
      if(pr) {
        zest <- rbind(time, surv)
        dimnames(zest) <- list(c("Time", "Survival"),
                               rep("", length(time)))
        if(slevp)cat("\nEstimates for ", slev[i], "\n\n")
        print(zest, digits=3)
      }
      if(conf.int > 0) {
        if(logt) {
          p <- pl(p, time, lower, type='lower', col=col[i], slev=sleva[i])
          p <- pl(p, time, upper, type='upper', col=col[i], slev=sleva[i])
        }
        else {
          p <- pl(p, c(min(time), time), c(fun(1), lower),
                  col=col[i], slev=slev[i], type='lower')  # see survplot ?max(tim)?
          p <- pl(p, c(min(time), time), c(fun(1), upper),
                  col=col[i], slev=slev[i], type='upper')
        }
      }
    }

    if(funtype %in% c('identity', 'inverse') && ns == 2 && conf.int > 0) {
      z <- survdiffplotp(fit.orig, fun=fun, conf.int=conf.int,
                         convert=convert, xlim=xlim, pobj=p)
      g <- plotly::toRGB('gray')
      p <- pl(p, z$time, z$lower, type='diff lower', col=g, slev='')
      p <- pl(p, z$time, z$upper, type='diff upper', col=g, slev='')
    }      

    
    if(aehaz || length(times)) {
      un <- if(units == ' ' | units == '') ''
            else
              paste('/', tolower(units), sep='')
      haz <- round(nevents / totaltime, 4)
      txt <- paste(nevents, 'events')
      if(aehaz) txt <- paste0(txt, '<br>', nbsp, nbsp, nbsp,
                              htmlGreek('lambda'), '=', 
                              haz, un, sep='')
      z <- paste(paste0(sleva, ':', txt), collapse='<br>')
      if(length(times)) {
        z2 <- paste0('<br>Cumulative<br>Incidence at<br>',
                     't=', paste(times, collapse=', '), ' ', units,
                     if(un !='') 's', '<br>',
                     paste0(sleva, ':', cuminc, collapse='<br>'))
        z <- paste0(z, z2)
        }
    }
    else z <- paste(paste0(sleva, ': ', nevents, ' events'), collapse='<br>')
    
    ## Add empty trace just to add to bottom of legend.  Used to have x=~NA y=~NA
    ## but plotly update made that point ignored in every way

    dam <- data.frame(x=xlim[1], y=ylim[1])
    p <- plotly::add_markers(p, x=~ x, y=~ y,               # mode='markers',
                           marker=list(symbol='asterisk'),  # suppresses pt
                           name=z, data=dam)
        
    xaxis <- list(range=xlim, title=xlab)
    if(! logt) xaxis <-
                 c(xaxis,
                   list(tickvals = seq(xlim[1], max(pretty(xlim)), time.inc)))
    
    plotly::layout(p,
                   xaxis=xaxis, 
                   yaxis=list(range=ylim, title=ylab))
}


