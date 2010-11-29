predict.lrm <- function(object, ..., 
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		"terms", "cterms", "ccterms", "adjto", "adjto.data.frame",
          "model.frame"),
		se.fit=FALSE, codes=FALSE) {

type <- match.arg(type)
if(type %nin% c("fitted","fitted.ind","mean"))
  return(predictrms(object,...,type=type, se.fit=se.fit))

xb <- predictrms(object, ..., type="lp", se.fit=FALSE)
rnam <- names(xb)
ns <- object$non.slopes
cnam <- names(object$coef[1:ns])
if(se.fit)warning('se.fit not supported with type="fitted" or type="mean"')
if(ns==1 & type=="mean")
  stop('type="mean" makes no sense with a binary response')
if(ns==1) return(1/(1+exp(-xb)))
intcept <- object$coef[1:ns]
xb <- xb - intcept[1]
xb <- sapply(intcept, "+", xb)
P <- 1/(1+exp(-xb))
nam <- names(object$freq)
if(is.matrix(P))dimnames(P) <- list(rnam, cnam)
else names(P) <- names(object$coef[1:ns])
if(type=="fitted") return(P)

#type="mean" or "fitted.ind"
vals <- names(object$freq)
P <- matrix(P,ncol=ns)
Peq <- cbind(1,P)-cbind(P,0)
if(type=="fitted.ind")
  {
    ynam <- as.character(attr(object$terms,"formula")[2])
    ynam <- paste(ynam,"=", vals, sep="")
    dimnames(Peq) <- list(rnam, ynam)
    return(drop(Peq))
  }

#type="mean"
if(codes) vals <- 1:length(object$freq)
else
  {
    vals <- as.numeric(vals)
    if(any(is.na(vals)))
      stop('values of response levels must be numeric for type="mean" and codes=F')
  }
m <- drop(Peq %*% vals)
names(m) <- rnam
m

}

Mean.lrm <- function(object, codes=FALSE, ...)
  {
    ns <- object$non.slopes
    if(ns < 2) stop('using this function only makes sense for >2 ordered response categories')
    if(codes) vals <- 1:length(object$freq)
    else
    {
      vals <- as.numeric(names(object$freq))
      if(any(is.na(vals)))
        stop('values of response levels must be numeric for codes=FALSE')
    }
    f <- function(lp=numeric(0), intercepts=numeric(0), values=numeric(0))
      {
        ns <- length(intercepts)
        lp <- lp - intercepts[1]
        xb <- sapply(intercepts, '+', lp)
        P  <- matrix(1/(1+exp(-xb)), ncol=ns)
        P  <- cbind(1,P) - cbind(P,0)
        m  <- drop(P %*% values)
        names(m) <- names(lp)
        m
      }
    formals(f) <- list(lp=numeric(0), intercepts=object$coef[1:ns],
                       values=vals)
    f
  }
