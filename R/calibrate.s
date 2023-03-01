calibrate <- function(fit, ...) UseMethod("calibrate")

print.calibrate <- function(x, B=Inf, ...)
{
  at <- attributes(x)
  predicted <- at$predicted
  dput(at$call)
  cat('\n\nn=', length(predicted), '  B=', at$B,
      '  u=', at$u, ' ', at$units, '\n\n', sep='')
  
  stratified <- 'KM' %in% colnames(x)
  if(stratified){
    attributes(x) <- at[c('dim','dimnames')]
    print.default(x)
  }
  else
    if(length(predicted)) {
      s <- !is.na(x[,'pred'] + x[,'calibrated.corrected'])
      err <- predicted -
        approxExtrap(x[s,'pred'],x[s,'calibrated.corrected'], 
                     xout=predicted, ties=mean)$y
      cat('\nMean |error|:', format(mean(abs(err))),
          '  0.9 Quantile of |error|:',
          format(quantile(err, 0.9, na.rm=TRUE)),
          '\n', sep='')
    }
  kept <- at$kept
  if(length(kept)) {
    cat("\nFactors Retained in Backwards Elimination\n\n")
    varin <- ifelse(kept, '*', ' ')
    print(varin[1:min(nrow(varin), B),], quote=FALSE)
    cat("\nFrequencies of Numbers of Factors Retained\n\n")
    nkept <- apply(kept, 1, sum)
    tkept <- table(nkept)
    names(dimnames(tkept)) <- NULL
    print(tkept)
  }
  invisible()
}

plot.calibrate <-
  function(x, xlab, ylab, subtitles=TRUE,
           conf.int=TRUE, cex.subtitles=.75,
           riskdist=TRUE, add=FALSE,
           scat1d.opts=list(nhistSpike=200),
           par.corrected=NULL,
           ...)
{
  at    <- attributes(x)
  u     <- at$u
  units <- at$units
  if(length(par.corrected) && ! is.list(par.corrected))
    stop('par.corrected must be a list')
  z <- list(col='blue', lty=1, lwd=1, pch=4)
  if(! length(par.corrected)) par.corrected <- z
  else for(n in setdiff(names(z), names(par.corrected)))
         par.corrected[[n]] <- z[[n]]
  
  predicted <- at$predicted
  if('KM' %in% colnames(x)) {
    type  <- 'stratified'
    pred  <- x[,"mean.predicted"]
    cal   <- x[,"KM"]
    cal.corrected <- x[,"KM.corrected"]
    se <- x[,"std.err"]
  }
  else {
    type <- 'smooth'
    pred <- x[,'pred']
    cal  <- x[,'calibrated']
    cal.corrected <- x[,'calibrated.corrected']
    se <- NULL
  }
  
  un <- if(u==1) paste(units, 's', sep='') else units
  if(missing(xlab)) xlab <- paste("Predicted ",format(u),units,"Survival")
  if(missing(ylab)) ylab <- paste("Fraction Surviving ",format(u)," ",un,
                                  sep="")

  ##Remember that groupkm stored the se of the log survival

  if(length(se) && conf.int) {
    ## Compute confidence limits for survival based on -log survival,
    ## constraining to be in [0,1]; d = std.error of cum hazard * z value
    ciupper <- function(surv, d) ifelse(surv==0, 0, pmin(1, surv*exp(d)))
    cilower <- function(surv, d) ifelse(surv==0, 0, surv*exp(-d))

    errbar(pred, cal, cilower(cal, 1.959964*se), ciupper(cal, 1.959964*se),
           xlab=xlab, ylab=ylab, type="b", add=add, ...)
  }
  else
    if(add) lines(pred, cal, type=if(type=='smooth') 'l' else 'b')
    else
      plot(pred, cal, xlab=xlab, ylab=ylab,
           type=if(type=='smooth')'l' else "b", ...)
  
  err <- NULL
  if(riskdist && length(predicted)) {
    do.call('scat1d', c(list(x=predicted), scat1d.opts))
    if(type=='smooth') {
      s <- !is.na(pred + cal.corrected)
      err <- predicted -
        approxExtrap(pred[s], cal.corrected[s], xout=predicted, ties=mean)$y
    }
  }
  
  if(subtitles && !add) {
    if(type=='smooth') {
      Col <- par.corrected$col
      substring(Col, 1, 1) <- toupper(substring(Col, 1, 1))
      title(sub=sprintf('Black: observed  Gray: ideal\n%s : optimism corrected',
              Col),
            adj=0, cex.sub=cex.subtitles)
      w <- if(length(err))
             paste('B=', at$B, ' based on ', at$what,
                   '\nMean |error|=', round(mean(abs(err)), 3),
                   '  0.9 Quantile=',
                   round(quantile(abs(err), .9, na.rm=TRUE), 3),
                   sep='')
           else
             paste('B=', at$B, '\nBased on ', at$what, sep='')
      title(sub=w, adj=1, cex.sub=cex.subtitles)
    }
    else {
      title(sub=paste("n=", at$n, " d=", at$d, " p=", at$p,
              ", ", at$m, " subjects per group\nGray: ideal", sep=""),
            adj=0, cex.sub=cex.subtitles)
      title(sub=paste("X - resampling optimism added, B=", at$B,
              "\nBased on ", at$what, sep=""),
            adj=1, cex.sub=cex.subtitles)
    }
  }
  
  abline(0, 1, col=gray(.9))	#ideal line
  if(type=='stratified')
    points(pred, cal.corrected, pch=par.corrected$pch, col=par.corrected$col)
  else
    lines (pred, cal.corrected, col=par.corrected$col, lty=par.corrected$lty,
           lwd=par.corrected$lwd)
  
  invisible()
}
