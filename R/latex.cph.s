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
  if(!length(which) & !inline)
    {
      if(length(strata)==0)
        {
          w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\} = S_{0}(t)^{{\\textstyle e}^{X\\beta}}, {\\rm \\ \\ where} \\\\ \\]",sep=""))
        }
      else
        {
          sname <- atr$name[atr$assume.code==8]
          strata.sub <- letters[8+(1:length(sname))]
          s <- paste("{\\rm ",sname,"}=",strata.sub,sep="")
          s <- paste(s, collapse=",")
          w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\ |\\ ",s,"\\}=S_{",
                         paste(strata.sub,collapse=""),
                         "}(t)^{{\\textstyle e}^{X\\beta}}, {\\rm \\ \\ where} \\\\ \\]", sep=""))
        }
    }
  if(!length(which)) which <- 1:length(atr$name)
  if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]
  if(! md) cat(w, sep=if(length(w))"\n" else "", file=file, append=append)

  Z <- latexrms(f, file=file, append=TRUE, which=which, varnames=varnames, 
                columns=columns, 
                before=before, after=after,
                prefix=if(!whichThere)"X\\hat{\\beta}" else NULL, 
                intercept=Intercept, inline=inline,
                pretrans=pretrans, digits=digits, size=size)
  if(md) Z <- c(paste0(w, '\n'), as.character(z))

  if(inline)
    return(if(md) htmltools::HTML(Z) else Z)
  
  ss <- f$surv.summary
  if(surv && length(ss)) {
    fs <- levels(f$strata)   # was f$strata
    nstrat <- 0; if(length(fs)) nstrat <- length(fs)
    times <- as.numeric(dimnames(ss)[[1]])
    maxtime <- f$maxtime
    if(max(times)>=maxtime) maxt <- FALSE
    if(nstrat==0) {
      s <- matrix(ss[, , 1], ncol=1)
      if(maxt) {
        s <- cbind(s, f$surv[L <- length(f$surv)])
        times <- c(times, f$time[L]) 
      }
      dimnames(s) <- list(format(times), "$S_{0}(t)$")
      if(md) {
        z <- htmlTable::txtRound(s, digits=dec)
        z <- htmlTable::htmlTable(z, rowlabel='$t$', escape.html=FALSE,
                                  css.cell='min-width: 9em;')
        Z <- c(Z, as.character(z))
      }
      else
        latex(s, file=file, append=TRUE, rowlabel="$t$",
              rowlabel.just="r",
              dec=dec, table.env=FALSE)
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

      dimnames(s) <- list(format(times),
                          paste("$S_{", n, "}(t)$", sep=""))
      if(md) {
        z <- htmlTable::txtRound(s, digits=dec)
        Z <- c(Z, as.character(
                    htmlTable::htmlTable(z, rowlabel='$t$',
                                         escape.html=FALSE,
                                         css.cell='min-width: 9em;')))
      }
      else
        latex(s, file=file, append=TRUE,
              rowlabel="$t$", rowlabel.just="r",
              dec=dec, table.env=FALSE)
    }
  }
  if(md) htmltools::HTML(Z)
}
