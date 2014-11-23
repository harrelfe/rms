datadist <- function(..., data, q.display, q.effect=c(.25,.75),
                     adjto.cat=c('mode','first'), n.unique=10)
{
  adjto.cat <- match.arg(adjto.cat)
  X <- list(...)

  argnames <- as.character(sys.call())[-1]
  
  if(inherits(x <- X[[1]],"datadist")) {
    Limits <- x$limits
    Values <- x$values
    X[[1]] <- NULL
    argnames <- argnames[-1]
  }
  else {
    Limits <- list()
    Values <- list()
  }
  
  if(is.data.frame(X[[1]])) {
    if(length(X) > 1) stop('when the first argument is a data frame, no other variables may be specified')
    X <- X[[1]]
  }
  
  else
    if(is.recursive(X[[1]]) &&
       length(Terms <- X[[1]]$terms) && length(D <- attr(Terms,"Design"))) {
      n <- D$name[D$assume != "interaction"]
      X <- list()
      if(missing(data))
        for(nm in n) X[[nm]] <- eval.parent(nm)
      else
        if(length(names(data))) {
          j <- match(n, names(data), 0)
          if(any(j == 0)) stop(paste("variable(s)",
                   paste(n[j == 0],collapse=" "),
                   "in model not found on data=, \nwhich has variables",
                   paste(names(data),collapse=" ")))
          for(nm in n) X[[nm]] <- data[[nm]]
        }
        else for(nm in n) X[[nm]] <- get(nm, data)  
    }
    else {
      if(length(X) & !length(names(X))) names(X) <- argnames[1 : length(X)]
      
### NEED TO FIX: R has no database.object
      if(!missing(data)) {
        ## This duplicative code is for efficiency for large data frames
        stop('program logic error')
            if(length(X)) {
              ## if(is.numeric(data)) X <- c(X,database.object(data))
              ## else
              X <- c(X, data)
            }
            else {
              ## if(is.numeric(data)) X <- database.object(data)
              ## else
              X <- data
            }
      }
    }
  nam <- names(X)
  p <- length(nam)
  if(p == 0) stop("you must specify individual variables or a data frame")
  
  maxl <- 0
  for(i in 1 : p) {
    values <- NULL
    x <- X[[i]]
    if(is.character(x)) x <- as.factor(x)
    lx <- length(x)
    lev <- levels(x)
    ll <- length(lev)
    limits <- rep(NA, 5)
    if(is.matrix(x) | (i > 1 && lx != maxl))
      warning(paste(nam[i],"is a matrix or has incorrect length; ignored"))
      else {
        if(ll && (ll < length(x))) values <- lev   # if # levels=length(x) is ID variable
        ## First look for ordered variable with numeric levels (scored() var)
        if(is.ordered(x) && all.is.numeric(lev)) {
          levx <- sort(as.numeric(lev))
          limits <- c(levx[1],levx[(ll+1)/2],levx[ll],levx[1],levx[ll],
                      levx[1],levx[ll])
          values <- levx
        }
        
        else if(ll) {
          adjto <- if(adjto.cat == 'first') lev[1] else {
            tab <- table(x)
            (names(tab)[tab == max(tab)])[1]
          }
          limits <- factor(c(NA,adjto,NA,lev[1],lev[ll],lev[1],lev[ll]),
                           levels=lev) 
          ## non-ordered categorical
        }
            else {	 	# regular numeric variable
              clx <- setdiff(class(x), c('integer', 'numeric'))
              ## Above prevents rounding of quantiles to integers
              y <- x[!is.na(x)]
              n <- length(y)
              if(n < 2)
                stop(paste("fewer than 2 non-missing observations for",nam[i]))
              values <- sort(unique(y))
              names(values) <- NULL
              nunique <- length(values)
              if(nunique < 2) {
                warning(paste(nam[i],"is constant"))
                limits <- rep(y[1], 7)
              }
              else {
                r <- range(values)
                limits[6 : 7] <- r
                if(nunique<4) q <- r else {
                  if(missing(q.display)) {
                    q.display <- 10 / max(n, 200)
                    q.display <- c(q.display, 1 - q.display)
                  }
                  q <- quantile(unclass(y), q.display)	}  #chron obj. not work here
                limits[4] <- q[1]; limits[5] <- q[2]
                ## check for very poorly distributed categorical numeric variable
                if(limits[4] == limits[5]) limits[4 : 5] <- r
                    
                ## Use low category if binary var, middle if 3-level, median otherwise
                if(nunique < 3) limits[2] <- values[1] else
                if(nunique == 3) limits[2] <- values[2] else
                limits[2] <- median(unclass(y))

                if(nunique < 4) q <- r else
                q <- quantile(unclass(y), q.effect)
                limits[1] <- q[1]; limits[3] <- q[2]
                if(limits[1] == limits[3]) limits[c(1,3)] <- r
                if(nunique > n.unique) values <- NULL
                class(limits) <- clx
              }
            }
        Limits[[nam[i]]] <- limits
        if(length(values)) Values[[nam[i]]] <- values
        maxl <- max(maxl, lx)
      }
  }
  
  Limits <- structure(Limits, class="data.frame", 
                      row.names=c("Low:effect","Adjust to",
                        "High:effect","Low:prediction",
                        "High:prediction","Low","High"))
  ##data.frame(Limits) gives error with chron objects
  
  d <- list(limits=Limits, values=Values)
  class(d) <- "datadist"
  d
}

print.datadist <- function(x, ...)
{
  lim <- x$limits
  for(n in names(lim)) {
    z <- lim[[n]]
    if(inherits(z,"dates") | inherits(z,"times"))
      lim[[n]] <- factor(format(z))
  }
  if(length(lim)) print(lim)
  ##print.data.frame doesn't print chron objects correctly
  if(length(V <- x$values)) {
    cat("\nValues:\n\n")
    wid <- .Options$width
    for(n in names(V)) {
      v <- V[[n]]
      if(length(v) == 0) next  # for gendata
      if(is.character(v) && length(v) > 80) 
        v <- c(v[1 : 20], paste("+", length(v), "others"))
      w <- if(is.character(v)) v else format(v)
      nc <- nchar(paste(w, collapse=" "))
      if(nc+nchar(n) + 4 > wid) {cat(n,":\n"); print(v, quote=FALSE)}
      else cat(n,":",w,"\n")
    }
  }
  invisible()
}
