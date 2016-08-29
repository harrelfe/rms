gIndex <-
  function(object, partials=TRUE, type=c('ccterms', 'cterms', 'terms'),
           lplabel=if(length(object$scale) && is.character(object$scale))
           object$scale[1] else 'X*Beta',
           fun, funlabel=if(missing(fun)) character(0) else
           deparse(substitute(fun)),
           postfun=if(length(object$scale)==2) exp else NULL,
           postlabel=if(length(postfun))
           ifelse(missing(postfun),
                  if((length(object$scale) > 1) &&
                     is.character(object$scale)) object$scale[2] else
                     'Anti-log',
                     deparse(substitute(postfun))) else character(0),
           ...)
{	
  obj.name <- as.character(sys.call())[2]
  type <- match.arg(type)
  labels <- attr(object, 'Design')$label

  lp <- predict(object, ...)

  if(partials)
    {
      terms <- predict(object, type=type)
      if(nrow(terms) != length(lp))
        warning('expected predicted linear predictors and terms to have same no. of rows')
      p <- ncol(terms)
      g <- matrix(0, nrow=p, ncol=1 + (length(postfun) > 0),
                  dimnames=list(colnames(terms),
                    c(lplabel, postlabel)))
      for(i in 1:p)
        {
          gmd   <- GiniMd(terms[,i], na.rm=TRUE)
          g[i,] <- c(gmd, if(length(postfun)) postfun(gmd)) 
        }
    }

  gmd <- GiniMd(lp, na.rm=TRUE)
  Total <- matrix(c(gmd, if(length(postfun)) postfun(gmd)),
                  nrow=1, ncol=1 + (length(postfun) > 0),
                  dimnames=list('Total', c(lplabel, postlabel)))
  g <- if(partials) rbind(g, Total) else Total
  gtrans <- NULL
  if(!missing(fun))
    {
      gtrans <- GiniMd(fun(lp), na.rm=TRUE)
      names(gtrans) <- funlabel
    }
  structure(g, gtrans=gtrans, class='gIndex',
            lplabel=lplabel, funlabel=funlabel,
            postlabel=postlabel, partials=partials,
            labels=c(labels, Total='Total'), type=type, formula=formula(object))
}

print.gIndex <-
  function(x, digits=4, abbrev=FALSE, vnames=c("names","labels"), ...)
{
  vnames <- match.arg(vnames)
  at <- attributes(x)
  
  if(vnames == 'labels')
    {
      lab <- at$labels[rownames(x)]
      rownames(x) <- if(abbreviate) abbreviate(lab) else lab
    }
  cat('\ng Index: ', format(at$formula), '\n\n')
  x <- matrix(x, nrow=nrow(x),
              dimnames=list(rownames(x), c(at$lplabel, at$postlabel)))
  print(x, digits=digits)
  if(length(gtrans <- at$gtrans))
    cat('\ng Index on transformed linear predictors (', names(gtrans), '): ',
        format(gtrans, digits=digits), '\n', sep='')
  cat('\n')
  invisible()
}

plot.gIndex <- function(x, what=c('pre', 'post'),
                          xlab=NULL, pch=16, rm.totals=FALSE,
                          sort=c('descending', 'ascending', 'none'), ...)
  {
      what <- match.arg(what)
      sort <- match.arg(sort)
      at <- attributes(x)

      if(!length(xlab))
        xlab <- paste('g Index:',
                      if(what=='pre') at$lplabel else at$postlabel)
      x <- if(what=='pre') x[, 1] else x[, 2]
      if(rm.totals) x <- x[-length(x)]

      x <- switch(sort,
               descending=-sort(-x),
               ascending=sort(x),
               none=x)

      dotchart3(x, xlab=xlab, pch=pch, ...)
      invisible(x)
  }
                          

