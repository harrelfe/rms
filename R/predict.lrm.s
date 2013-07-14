predict.lrm <- function(object, ..., 
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		"terms", "cterms", "ccterms", "adjto", "adjto.data.frame",
          "model.frame"),
		se.fit=FALSE, codes=FALSE)
{
  type <- match.arg(type)
  if(type %nin% c("fitted","fitted.ind", "mean"))
    return(predictrms(object,...,type=type, se.fit=se.fit))

  xb <- predictrms(object, ..., type="lp", se.fit=FALSE)
  rnam <- names(xb)
  ns <- object$non.slopes
  cnam <- names(object$coef[1:ns])
  trans <- object$trans
  ## If orm object get cumulative probability function used
  cumprob <- if(length(trans)) trans$cumprob else plogis
  if(se.fit)
    warning('se.fit not supported with type="fitted" or type="mean"')
  if(ns == 1 & type == "mean")
    stop('type="mean" makes no sense with a binary response')
  if(ns == 1) return(cumprob(xb))
  intcept <- object$coef[1:ns]
  interceptRef <- object$interceptRef
  if(!length(interceptRef)) interceptRef <- 1
  xb <- xb - intcept[interceptRef]
  xb <- sapply(intcept, "+", xb)
  P <- cumprob(xb)
  nam <- names(object$freq)
  if(is.matrix(P)) dimnames(P) <- list(rnam, cnam)
  else names(P) <- names(object$coef[1:ns])
  if(type=="fitted") return(P)

  ##type="mean" or "fitted.ind"
  vals <- names(object$freq)
  P   <- matrix(P, ncol=ns)
  Peq <- cbind(1, P) - cbind(P, 0)
  if(type == "fitted.ind") {
    ynam <- as.character(attr(object$terms, "formula")[2])
    ynam <- paste(ynam, "=", vals, sep="")
    dimnames(Peq) <- list(rnam, ynam)
    return(drop(Peq))
  }
  
  ##type="mean"
  if(codes) vals <- 1:length(object$freq)
  else {
    vals <- as.numeric(vals)
    if(any(is.na(vals)))
      stop('values of response levels must be numeric for type="mean" and codes=F')
  }
  m <- drop(Peq %*% vals)
  names(m) <- rnam
  m
}

predict.orm <- function(object, ..., 
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		"terms", "cterms", "ccterms", "adjto", "adjto.data.frame",
          "model.frame"),
		se.fit=FALSE, codes=FALSE)
{
  type <- match.arg(type)
  predict.lrm(object, ..., type=type, se.fit=se.fit, codes=codes)
}

Mean.lrm <- function(object, codes=FALSE, ...)
  {
    ns <- object$non.slopes
    if(ns < 2)
stop('using this function only makes sense for >2 ordered response categories')
    if(codes) vals <- 1:length(object$freq)
    else {
      vals <- object$yunique
      if(!length(vals)) vals <- names(object$freq)
      vals <- as.numeric(vals)
      if(any(is.na(vals)))
        stop('values of response levels must be numeric for codes=FALSE')
    }
    f <- function(lp=numeric(0), intercepts=numeric(0), values=numeric(0),
                  interceptRef=integer(0), cumprob=cumprob)
      {
        ns <- length(intercepts)
        lp <- lp - intercepts[interceptRef]
        xb <- sapply(intercepts, '+', lp)
        P  <- matrix(cumprob(xb), ncol=ns)
        P  <- cbind(1, P) - cbind(P, 0)
        m  <- drop(P %*% values)
        names(m) <- names(lp)
        m
      }
    cumprob <- object$trans$cumprob
    if(!length(cumprob)) cumprob <- plogis
    ir <- object$interceptRef
    if(!length(ir)) ir <- 1
    formals(f) <- list(lp=numeric(0), intercepts=object$coef[1:ns],
                       values=vals, interceptRef=ir, cumprob=cumprob)
    f
  }
