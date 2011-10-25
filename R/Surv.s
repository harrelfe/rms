Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2'),
		       origin=0)
{
  nam  <- as.character(sys.call())[-1]
  mtype <- missing(type)
  type <- match.arg(type)
 
  
  ng  <- (!missing(time)) + (!missing(time2)) + (!missing(event))
  if (missing(type))
    {
      if (ng==1 || ng==2) type <- 'right'
      else if (ng==3)     type <- 'counting'
      else stop("Invalid number of arguments")
	}
  if(ng==1)
    {
      tvar <- time
      svar <- NULL
      nam  <- nam[1]
    }
  else if (type=='right' || type=='left')
    {
      tvar <- time
      svar <- time2
      nam <- nam[1:2]
    }
  else
    {
      tvar <- time2
      svar <- event
      nam  <- nam[2:3]
    }

  g <- list()
  if(!missing(time )) g$time  <- time
  if(!missing(time2)) g$time2 <- time2
  if(!missing(event)) g$event <- event
  if(!mtype)          g$type  <- type
  if(!missing(origin)) g$origin <- origin
  
  surv <- survival:::Surv
  ss <- do.call('surv', g)

  uni <- valueUnit(tvar)
  if(!length(uni)) uni <- "Day"
  tlab <- attr(tvar, "label")
  if(!length(tlab)) tlab <- nam[1]
  elab <- attr(svar, "label")
  if(!length(elab) && length(nam)>1) elab <- nam[2]
  valueUnit(ss) <- uni
  attr(ss,"time.label") <- tlab
  attr(ss,"event.label") <- elab
  ss
}

"[.Surv" <- function(x, ..., drop=FALSE)
{
  atr <- attributes(x)
  atr$dim <- NULL; atr$dimnames <- NULL
  if (missing(..2)) {
    cl <- class(x)
    class(x) <- NULL
    x <- NextMethod('[')
    class(x) <- cl
    attributes(x) <- c(attributes(x), atr)
	x
  }
  else {
    oldClass(x) <- NULL
    NextMethod("[")
  }
}

