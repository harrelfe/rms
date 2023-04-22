latex.cph <-
  function(object, title, 
           file='',
           append=FALSE, surv=TRUE, maxt=FALSE, which=NULL, varnames, 
           columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           dec=3, pretrans=TRUE,
           caption=NULL, digits=.Options$digits, size='', ...)
{
  md <- prType() %in% c('html', 'md', 'markdown')
  filegiven <- file != ''
  f <- object
  whichThere <- length(which)
  
  atr <- f$Design
  
  lev <- names(f$freq)
  Intercept <- -f$center
  strata <- levels(f$strata)    ## was f$strata
  w <- if(length(caption)) {
         if(md) paste('<div align=center><strong>', caption,
                      '</strong></div>')
         else
           paste('\\begin{center} \\bf',caption,'\\end{center}')
         }
  if(! length(which) & !inline) {
    if(length(strata)==0)
      w <- c(w, paste("$$\\Pr(T\\geq t~|~X) = S_{0}(t)^{\\mathrm{e}^{X\\beta}},~~ \\mathrm{where}$$",sep=""))
    else {
      sname <- atr$name[atr$assume.code==8]
      strata.sub <- letters[8 + (1 : length(sname))]
      s <- paste("\\mathrm{",sname,"}=",strata.sub,sep="")
      s <- paste(s, collapse=",")
      w <- c(w,paste("$$\\Pr(T\\geq t~|~X,",s,")=S_{",
                     paste(strata.sub,collapse=""),
                     "}(t)^{\\mathrm{e}^{X\\beta}},~~\\mathrm{where}$$", sep=""))
    }
  }
  if(!length(which)) which <- 1:length(atr$name)
  if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]
  ltx <-
    latexrms(f, file='', append=TRUE, which=which, varnames=varnames, 
             columns=columns, 
             before=before, after=after,
             prefix=if(!whichThere)"X\\hat{\\beta}" else NULL, 
             intercept=Intercept, inline=inline,
             pretrans=pretrans, digits=digits, size=size)
  if(inline) return(ltx)
  w <- c(w, ltx)

  htmlTab <- function(s) {
    s <- cbind('$t$'= as.numeric(rownames(s)), s)
    for(j in 1 : ncol(s)) s[, j] <- round(s[, j], dec)
    if (requireNamespace("kableExtra", quietly=TRUE)) {
      as.character(
        knitr::kable(s, format='html',
                     align='r', row.names=FALSE) |>
        kableExtra::kable_styling(full_width=FALSE)   )
    } else {
      as.character(
        knitr::kable(s, format='html',
                     align='r', row.names=FALSE) )
    }
  }
  
  ss <- f$surv.summary
  if(surv && length(ss)) {
    tf <- tempfile()
    fs <- levels(f$strata)
    nstrat <- 0; if(length(fs)) nstrat <- length(fs)
    times <- as.numeric(dimnames(ss)[[1]])
    maxtime <- f$maxtime
    if(max(times) >= maxtime) maxt <- FALSE
    if(nstrat == 0) {
      s <- matrix(ss[, , 1], ncol=1)
      if(maxt) {
        s <- cbind(s, f$surv[L <- length(f$surv)])
        times <- c(times, f$time[L]) 
      }
      dimnames(s) <- list(t=format(times), "$S_{0}(t)$")
      if(md) {
#        z <- htmlTable::txtRound(s, digits=dec)
#        z <- htmlTable::htmlTable(z, rowlabel='$t$', escape.html=FALSE,
#                                  css.cell='min-width: 9em;')
        w <- c(w, htmlTab(s))
      } else {
        latex(s, file=tf, append=TRUE, rowlabel="$t$",
              rowlabel.just="r",
              dec=dec, table.env=FALSE)
        w <- c(w, readLines(tf))
        }
    } else {
          
      ## Change . to ,blank
      n <- sedit(paste(fs,',',sep=''), '.', ', ')
      ## Change sname=*, to *,
      n <- sedit(n, paste(sname,'=*,',sep=''), rep('*, ', length(sname)))
      n <- substring(n, 1, nchar(n) - sum(atr$assume.code == 8) - 1)
      s <- ss[, , 1]
      if(maxt) {
        smax <- rep(NA, nstrat)
        for(i in 1 : nstrat)
          smax[i] <- f$surv[[i]][abs(f$time[[i]]-maxtime) < 0.001]
        s <- rbind(s, smax)
        times <- c(times, maxtime)
      }    

      dimnames(s) <- list(t=format(times),
                          paste("$S_{", n, "}(t)$", sep=""))
      if(md) {
#        z <- htmlTable::txtRound(s, digits=dec)
#        z <- htmlTable::htmlTable(z, rowlabel='$t$',
#                                  escape.html=FALSE,
#                                  css.cell='min-width: 9em;')
        w <- c(w, htmlTab(s))
      }
      else {
        ltx < latex(s, file=tf, append=TRUE,
                    rowlabel="$t$", rowlabel.just="r",
                    dec=dec, table.env=FALSE)
        w <- c(w, readLines(tf))
      }
    }
  }
  if(filegiven || prType() == 'plain') {
    cat('\n', w, sep='\n', file=file, append=append)
    return(invisible())
  }
  rendHTML(w, html=FALSE)
}
