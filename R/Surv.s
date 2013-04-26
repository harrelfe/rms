Srv <- function(time, time2, event,
                type=c('right', 'left', 'interval', 'counting', 'interval2',
                  'mstate'), origin=0)
{
  nam   <- as.character(sys.call())[-1]
  mtype <- missing(type)
  type  <- match.arg(type)
   
  ng  <- (!missing(time)) + (!missing(time2)) + (!missing(event))
  if (mtype || type == 'mstate') {
    if (ng==1 || ng==2) type <- 'right'
    else if (ng==3)     type <- 'counting'
    else stop("No time variable!")
	}
  if(ng==1) {
    tvar <- time
    svar <- NULL
    tnam  <- deparse(substitute(time))
    snam <- ''
  }
  else if (type=='right' || type=='left') {
    tvar <- time
    svar <- if(missing(event)) time2 else event
    tnam <- deparse(substitute(time))
    snam <- if(missing(event)) deparse(substitute(time2)) else
                               deparse(substitute(event))
  }
  else {
    tvar <- time2
    svar <- if(!missing(event)) event
    tnam <- deparse(substitute(time2))
    snam <- if(!missing(event)) deparse(substitute(event)) 
  }

  g <- list()
  if(!missing(time ))  g$time   <- time
  if(!missing(time2))  g$time2  <- time2
  if(!missing(event))  g$event  <- event
  if(!mtype)           g$type   <- type
  if(!missing(origin)) g$origin <- origin
  
  ss <- do.call('Surv', g)

  uni <- valueUnit(tvar)
  if(!length(uni)) uni <- "Day"
  tlab <- attr(tvar, "label")
  if(!length(tlab)) tlab <- tnam
  elab <- attr(svar, "label")
  if(!length(elab)) elab <- snam
  valueUnit(ss) <- uni
  if(length(tlab) && tlab != '') attr(ss,"time.label")  <- tlab
  if(length(elab) && elab != '') attr(ss,"event.label") <- elab
  class(ss) <- c('Srv', class(ss))
  ss
}

"[.Srv" <- function(x, ..., drop=FALSE)
{
  atr <- attributes(x)
  atr$dim <- NULL; atr$dimnames <- NULL
  if (missing(..2)) {
    cl <- class(x)
    class(x) <- NULL
    x <- NextMethod('[')
    attributes(x) <- c(attributes(x), atr)
    class(x) <- cl
	x
  }
  else {
    class(x) <- NULL
    NextMethod("[")
  }
}

