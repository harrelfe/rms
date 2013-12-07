# design.trans   FEH 4 Oct 90
# Contains individual functions for creating sub-design matrices from
# vectors, for use with design().
# code  name
# 1	asis	leave variable coded as-is, get default name, label,
#		limits, values
# 2	pol	polynomial expansion
# 3	lsp	linear spline
# 4	rcs	restricted cubic spline
# 5	catg	category
# 7	scored	scored ordinal variable
# 8	strat	stratification factor
#10	matrx	matrix factor - used to keep groups of variables together
#		as one factor
#
#	des.args generic function for retrieving arguments
#	set.atr generic function to set attributes of sub design matrix
#	options sets default options
#	[.rms subsets variables, keeping attributes
#	gparms	retrieve parms for design or fit object.  Not used by any
#		of these routines, but used by analyst to force a new
#		fit to use same parms as previous fit for a given factor.
#	value.chk
#		Check a given list of values for a factor for validity,
#		or if list is NA, return list of possible values
#	
# Default label is attr(x,"label") or argument name if label= omitted.
# First argument can be as follows, using asis as an example:
#	asis(x, ...)		name="x", label=attr(x,"label") or "x"
#				if NULL
#	asis(w=abs(q), ...)	name="w", label=attr(x,"label") or "w"
#	asis(age=xx)		name="age", label=label attr or "age"
#	asis(x,label="Age, yr")	name="x", label="Age, yr"
#	asis(age=x,label=	name="age", label="Age in Years"
#		"Age in Years")
#	matrx(dx=cbind(dx1=dx1,dx2=dx2))	name="dx", individual names
#						dx1 and dx2
# For matrx, default label is list of column names.
# An additional argument, name, can be used to instead specify the name of the
# variable.  This is used when the functions are implicitly called from within
# design().
#
# The routines define dimnames for the returned object with column
# names = expanded list of names based on original name.
# assume.code is added to attributes of returned matrix.  Is 1-8
# corresponding to transformation routines asis-strat above, 10 for matrx.
# Adds attribute nonlinear, one element/column of expanded design matrix.
# nonlinear=T if column is a nonlinear expansion of original variable,
# F if linear part or not applicable
# (e.g. dummy variable for category -> F).  For matrx, all are linear.
#
# System options used: nknots for default number of knots in restr. cubic spline
# and poly.degree, default degree of polynomials 
# Second argument to routines is the parameters (parms) of the
# transformation (except for asis), defined as follows:
#
#	poly	order of polynomial, e.g. 2 for quadratic
#	lsp	list of knots
#	rcs	number of knots if parms=1 element (-> compute default
#		knot locations), actual knot locations if >2 elements
#		(2 knots not allowed for restr. cubic spline)
#	catg	list of value labels corresponding to values 1,2,3,...
#	scored	list of unique values of the variable
#	strat	list of value labels corresponding to values 1,2,3
#
# For catg and strat, parms are omitted if the variable is character or
# is already an S category variable.
#
# Argument retrieval: After variable and optional parms, other variables
# may be named or positional, in the following order: label, name.
# For matrx, parms are not allowed.
#
# Function to return list with elements name, parms, label. 
# corresponding to arguments in call to asis, etc.  parms=NULL if 
# parms.allowed=F.  Reason for going to this trouble is that first arg to
# asis, etc. is allowed to be a named argument to set a new name for it.
# With ordinary argument fetching, remaining arguments would have to be
# named.  This logic allows them to be named or positional in order:
# parms (if allowed), label.
#
# If options(Design.attr) is non-null, looks up attributes in elements
# in Design.attr corresponding to the name of the current variable.
# This is used to get predicted values when the original fitting
# function (e.g., rcs) derived parms of the transformation from the data.
#
des.args <- function(x,parms.allowed,call.args) {
  nam <- names(x)
  if(!length(nam)) nam <- rep("",5)
  name <- nam[1]
  if(name=="") {
    form <- formula(call("~",as.name("...y..."),call.args[[2]]))
    name <- var.inner(form)
  }
  pa <- parms.allowed
  argu <- function(x,karg, arg.name, parms.all, nm)	{
	if(!parms.all) karg <- karg-1
	k <- charmatch(arg.name,nm,0)	#k>0 : named arg found
    ## Added karg <= length(x) 9Apr02 for R; R doesn't return NULL
    ## like S+
	if(k>0) x[[k]] else 
	if(length(nm) < karg || nm[karg]!="") NULL else
     if(karg <= length(x)) x[[karg]] else NULL
  }
  if(parms.allowed) parms <- argu(x,2,"parms",pa,nam) else {
	parms <- NULL
	if(charmatch("parms",nam,0)>0)
      stop(paste("parms not allowed for",as.character(call.args[1])))
  }
 
  nm <- argu(x,5,"name",pa,nam)
  if(!is.null(nm)) name <- nm
  if(!is.null(.Options$Design.attr)) {
	atr <- .Options$Design.attr
	i <- charmatch(name, atr$name, 0)
	if(is.null(i))stop("program logic error for options(factor.number)")
	parmi <- atr$parms[[name]]
	return(list(name=atr$name[i],parms=parmi,label=atr$label[i],
                units=atr$units[i]))		# added units 9Jun99
  }

  label <- argu(x,3,"label",pa,nam)
  atx <- attributes(x[[1]])  # 9Jun99
  if(is.null(label)) label <- atx$label   # 9Jun99 attr(x[[1]],"label")
  if(is.null(label)) label <- name

  list(name=name,parms=parms,label=label,units=atx$units)  #9Jun99
  
}

## Function to list all attributes of new sub-design matrix
set.atr <- function(xd, x, z, colnames, assume, code, parms, nonlinear) {
  ##Note: x argument isn't used
  if(is.matrix(xd))
    list(dim=dim(xd),dimnames=list(NULL,colnames),class="rms",
         name=z$name, label=z$label, assume=assume, assume.code=code,
         parms=parms, 
         nonlinear=nonlinear,colnames=colnames,units=z$units)
  else list(dim=dim(xd), class="rms",
            name=z$name, label=z$label, assume=assume, assume.code=code,
            parms=parms, 
            nonlinear=nonlinear,colnames=colnames,units=z$units)
}

## asis transformation - no transformation	
asis <- function(...) {
  cal <- sys.call()
  xx <- list(...)
  z <- des.args(xx, FALSE, cal)
  xd <- xx[[1]]
  if(is.factor(xd)) {
    attr(xd,"class") <- NULL
  }
  if(!(is.numeric(xd) | is.logical(xd))) {
    stop(paste(z$name,"is not numeric"))
  }
  
  attributes(xd) <- set.atr(xd,xd,z,z$name,"asis",1,NULL,FALSE)
  xd
}

## matrx transformation - no transformation, keep original vars as matrix
## column names as parameter names, parms=column medians
matrx <- function(...) {

  cal <- sys.call()
  xx <- list(...)
#  prn(xx);prn(cal)
  z <- des.args(xx, FALSE, cal)
#  prn(z)

  xd <- xx[[1]]
  nc <- ncol(xd)
  if(!is.matrix(xd)) {
    stop(paste(z$name, "is not a matrix"))
  }
  colname <- dimnames(xd)[[2]]
  if(length(colname)==0 && nc > 0)
    colname <- paste(z$name, '[', 1:nc, ']', sep="")
  else if(z$label==z$name)
    z$label <- paste(colname, collapse=",")

  parms <- rep(NA, max(1, nc))
  if(length(xd)) for(i in 1:nc)
    parms[i] <- median(xd[,i], na.rm=TRUE)

  attributes(xd) <- set.atr(xd, NULL, z, colname, "matrix", 10, parms,
                            rep(FALSE,nc))
  xd
}

## Polynomial expansion
pol <- function(...) {

  cal <- sys.call()
  xx <- list(...)
  z <- des.args(xx,TRUE,cal)
  x <- xx[[1]]
  if(!is.numeric(x)) {
    stop(paste(z$name,"is not numeric"))
  }
  poly.degree <- .Options$poly.degree
  if(is.null(poly.degree)) {
    poly.degree <- 2
  }

  if(!is.null(z$parms)) {
    poly.degree <- z$parms
  }

  if(poly.degree<2){
    stop("order for polynomial must be 2,3,...")
  }

  xd <- matrix(1,nrow=length(x),ncol=poly.degree)
  nam <- z$name
  name <- character(poly.degree)
  name[1] <- nam
  xd[,1] <- x

  for(j in 2:poly.degree) {
    name[j] <- paste(nam,"^",j,sep="")
    xd[,j] <- x^j
  }

  attributes(xd) <- set.atr(xd,x,z,name,"polynomial",2,poly.degree,
                            c(FALSE,rep(TRUE,poly.degree-1)))
  xd
}


## Linear spline expansion
lsp <- function(...) {

  cal <- sys.call()
  xx <- list(...)
  z <- des.args(xx,TRUE,cal)
  x <- xx[[1]]
  if(!is.numeric(x)) {
    stop(paste(z$name,"is not numeric"))
  }
  parms <- z$parms
  if(is.null(parms) || any(is.na(parms)))  {
    stop("must specify knots for linear spline")
  }

  suffix <- NULL
  nam <- z$name
  lp <- length(parms)
  xd <- matrix(double(1),nrow=length(x),ncol=lp+1)
  name <- character(lp+1)
  xd[,1] <- x
  name[1] <- nam

  for(j in 1:lp) {
    suffix <- paste(suffix,"'",sep="")
    name[j+1] <- paste(nam,suffix,sep="")
    xd[,j+1] <- pmax(x-parms[j],0)
  }

  attributes(xd) <- set.atr(xd,x,z,name,"lspline",3,parms,c(FALSE,rep(TRUE,lp)))
  xd
}

## Restricted cubic spline expansion
rcs <- function(...) {
  
  cal <- sys.call()
  xx <- list(...)
  z <- des.args(xx, TRUE, cal)
  x <- xx[[1]]
  if(!is.numeric(x)) stop(paste(z$name, "is not numeric"))

  nknots <- .Options$nknots
  if(!length(nknots)) nknots <- 5

  parms <- z$parms
  if(!length(parms)) parms <- nknots

  if(length(parms)==1) {
    nknots <- parms
    knots <- NULL
    if(nknots == 0) {
      attributes(x) <- set.atr(x, x, z, z$name, "asis", 1, NULL, FALSE)
      return(x)
    }
  }
  else {
    nknots <- length(parms)
    knots <- parms
  }
  
  pc <- length(.Options$rcspc) && .Options$rcspc
  fractied <- .Options$fractied
  if(!length(fractied)) fractied <- 0.05
  
  if(!length(knots)) {
    xd <- rcspline.eval(x, nk=nknots, inclx=TRUE, pc=pc, fractied=fractied)
    knots <- attr(xd,"knots")
  }
  else xd <- rcspline.eval(x, knots=knots, inclx=TRUE, pc=pc, fractied=fractied)

  parms  <- knots
  nknots <- length(parms)
  nam    <- z$name
  primes <- paste(rep("'",nknots-1), collapse="")
  name   <- if(pc)
    paste(nam, substring(primes, 1, 1:(nknots-1)), sep="")
  else c(nam, paste(nam, substring(primes, 1, 1:(nknots-2)), sep=""))

  if(pc) attr(parms, 'pcparms') <- attr(xd, 'pcparms')
  attributes(xd) <-
    set.atr(xd, x, z, name, "rcspline", 4, parms,
            if(pc) rep(TRUE, nknots-1) else c(FALSE,rep(TRUE,nknots-2)))
  xd
}

## Category variable
catg <- function(...) {
  cal <- sys.call()
  xx <- list(...)
  z <- des.args(xx,TRUE,cal)
  nam <- z$name
  y <- xx[[1]]
  parms <- z$parms

  if(is.null(parms) & is.factor(y)) parms <- levels(y)

  if(is.null(parms)) {
    if(is.character(y)) {
      parms <- sort(unique(y[y != "" & y != " "]))
    } else {
      parms <- as.character(sort(unique(y[!is.na(y)])))
    }
  }
  
  if(!is.factor(y)) {
    x <- factor(y, levels=parms)
  } else {
    x <- y
  }

  if((is.character(y) && any(y!="" & y!=" " & is.na(x))) ||
     (is.numeric(y) & any(!is.na(y) & is.na(x)))) {
    stop(paste(nam,"has non-allowable values"))
  }

  if(all(is.na(x))) {
    stop(paste(nam,"has no non-missing observations"))
  }

  lp <- length(parms)
  if(lp < 2) stop(paste(nam,"has <2 category levels"))

  attributes(x) <- list(levels=parms,class=c("factor","rms"),
                        name=nam,label=z$label,assume="category",assume.code=5,
                        parms=parms,nonlinear=rep(FALSE,lp-1),
                        colnames=paste(nam,"=",parms[-1],sep=""))
  x
}

## Scored expansion   parms=unique values
scored <- function(...) {

  cal <- sys.call()
  xx <- list(...)
  z <- des.args(xx,TRUE,cal)
  parms <- z$parms
  nam <- z$name
  x <- xx[[1]]
  if(is.factor(x)) {
    levx <- as.numeric(levels(x))
    if(any(is.na(levx))) stop(paste("levels for",nam,"not numeric"))
    if(is.null(parms)) parms <- levx

    ## .Options$warn <- -1   #suppress warning about NAs
    oldopt <- options(warn=-1)
    on.exit(options(oldopt))
    x <- levx[x]
  }

  if(!is.numeric(x)) stop(paste(nam,"is not a numeric variable"))

  y <- sort(unique(x[!is.na(x)]))
  if(is.null(parms)) parms <- y

  parms <- sort(parms)
  n.unique <- length(parms)
  if(n.unique < 3) {
    stop("scored specified with <3 levels")
  }
  lp <- length(parms)-1

  ## Form contrast matrix of the form linear | dummy | dummy ...

  xd <- matrix(double(1), nrow=length(y), ncol=lp)
  xd[,1] <- y
  name <- character(lp)
  name[1] <- nam
  i <- 1
  for(k in parms[3:length(parms)]) {
    i <- i+1
    name[i] <- paste(nam,"=",k,sep="")
    xd[,i] <- y==k
  }

  dimnames(xd) <- list(NULL, name)

  x <- ordered(x)
  class(x) <- c("ordered","factor","rms")

  attributes(x) <- c(attributes(x),
                     list(name=nam,label=z$label,assume="scored",assume.code=7,
                          parms=parms,
                          nonlinear=c(FALSE,rep(TRUE,lp-1)), colnames=name,
                          contrasts=xd))
  x
}

## strat parms=value labels
strat <- function(...) {

  cal <- sys.call()
  xx <- list(...)
  y <- xx[[1]]
  z <- des.args(xx,TRUE,cal)
  parms <- z$parms
  if(is.null(parms)) parms <- levels(y)
  if(is.null(parms)) {
    if(is.character(y)) {
      parms <- sort(unique(y[y!="" & y!=" "]))
    } else parms <- as.character(sort(unique(y[!is.na(y)])))
  }

  nam <- z$name
  if(!is.factor(y)) {
    x <- factor(y,levels=parms)
  } else x <- y

  if((is.character(y) & any(y!="" & y!=" " & is.na(x))) ||
     (is.numeric(y) & any(!is.na(y) & is.na(x)))) {
    stop(paste(nam," has a non-allowable value"))
  }

  name <- nam
  attributes(x) <- list(levels=parms,class=c("factor","rms"),
                        name=nam, label=z$label, assume="strata", assume.code=8,
                        parms=parms, nonlinear=FALSE)
  x
}

## Function to subscript a variable, keeping attributes
## Is similar to [.smooth, but does not keep attribute NAs
"[.rms" <- function(x, ..., drop = FALSE) {
  ats <- attributes(x)
  ats$dimnames <- NULL
  ats$dim <- NULL
  ats$names <- NULL
  class(x) <- NULL
  y <- x[..., drop = drop]
  attributes(y) <- c(attributes(y), ats)
  y
}

## Function to get parms of factor in fit or design object "fit" with name
## given by second argument (without quotes)
gparms <- function(fit,...) {
  name <- as.character(sys.call())[3]
  atr <- fit$Design
  atr$parms[[name]]
}

## value.chk - if x=NA, returns list of possible values of factor i defined
##	in object f's attributes.  For continuous factors, returns n values
##	in default prediction range.  Use n=0 to return trio of effect
##	limits.  Use n<0 to return pretty(plotting range,nint=-n).
##       If type.range="full" uses the full range instead of default plot rng.
## If x is not NA, checks that list to see that each value is allowable
##	for the factor type, and returns x
## Last argument is object returned from Getlim (see Design.Misc)
## First argument is Design list

value.chk <- function(f, i, x, n, limval, type.range="plot")
{
  as     <- f$assume.code[i]
  name   <- f$name[i]
  parms  <- f$parms[[name]]
  isna   <- length(x)==1 && is.na(x)
  values <- limval$values[[name]]
  charval <- !is.null(values) && is.character(values)
  if(isna & as!=7) {
    if(is.null(limval) || match(name, dimnames(limval$limits)[[2]], 0)==0 ||
       is.na(limval$limits["Adjust to",name]))
      stop(paste("variable",name,"does not have limits defined by datadist"))
    
    limits <- limval$limits[,name]
    lim    <- if(type.range=="full") limits[6:7] else limits[4:5]
  }

  if(as<5 | as==6) {
      if(isna) {
        if(is.null(values)) {
          if(n==0) x <- limits[1:3]
          else {
            if(n>0) x <- seq(unclass(lim[1]), #handles chron
                             unclass(lim[2]),length=n)
            else x <- pretty(unclass(lim[1:2]), n=-n)
            class(x) <- class(lim)
          }
        } else x <- values
      } else {
        if(is.character(x) && !charval)
          stop(paste("character value not allowed for variable",
                     name))   #Allow any numeric value
        if(charval) {
          j <- match(x, values, 0)
          if(any(j==0))
            stop(paste("illegal values for categorical variable:",
                       paste(x[j==0],collapse=" "),"\nPossible levels:",
                       paste(values,collapse=" ")))
        }	
      }
    } else if(as==5|as==8) {
      if(isna) x <- parms
      else {
        j <- match(x, parms, 0)  #match converts x to char if needed
        if(any(j==0))
          stop(paste("illegal levels for categorical variable:",
                     paste(x[j==0],collapse=" "),"\nPossible levels:",
                     paste(parms,collapse=" ")))
        x
      }
    }
    else if(as==7) {
      if(isna) x <- parms
      else if(is.character(x))
        stop(paste("character value not allowed for",
                   "variable",name))
      else {
        j <- match(x, parms, 0)
        if(any(j==0)) {
          stop(paste("illegal levels for categorical variable:",
                     paste(x[j==0],collapse=" "),"\n","Possible levels:",
                     paste(parms,collapse=" ")))
        }
      }
    }
  
  invisible(x)
}

##ia.operator.s - restricted interaction operators for use with Design
##F. Harrell  8 Nov 91

##Set up proper attributes for a restricted interaction for a model
##such as y ~ rcs(x1) + rcs(x2) + x1 %ia% x2 or x1 %ia% rcs(x2)
##or rcs(x1) %ia% x2

"%ia%" <- function(x1, x2) {
  a1 <- attributes(x1)
  a2 <- attributes(x2)
  nam <- as.character(sys.call())[-1]

  redo <- function(x,nam) {
      if(is.null(attr(x,"assume.code"))) {
        if(!is.null(class(x)) && class(x)[1]=="ordered")
          x <- scored(x, name=nam)
        else if(is.character(x) | is.factor(x))
          x <- catg(x, name=nam)
        else if(is.matrix(x)) x <- matrx(x, name=nam)
        else x <- asis(x, name=nam)
      }
      ass <- attr(x,"assume.code")
      nam <- attr(x,"name")
      
      if(ass==5) {
        colnames <- attr(x,"colnames")
        len <- length(attr(x,"parms"))-1
      } else if(ass==8) {
        prm <- attr(x,"parms")
        colnames <- paste(nam,"=",prm[-1],sep="")
        len <- length(prm)-1
      } else if(ass==7) {
        prm <- attr(x,"parms")
        colnames <- c(nam,paste(nam,"=",prm[-(1:2)],sep=""))
        len <- length(prm)-1
      } else {
        if(is.null(ncol(x))) {
          len <- 1
          colnames <- nam
        } else {
          colnames <- dimnames(x)[[2]]
          len <- ncol(x)
        }
      }

      attr(x,"colnames") <- colnames
      attr(x,"len") <- len
      if(ass==8) attr(x,"nonlinear") <- rep(FALSE, len)
      x
    }

  x1 <- redo(x1,nam[1])
  x2 <- redo(x2,nam[2])
  a1 <- attributes(x1)
  a2 <- attributes(x2)
  n1 <- a1$colnames
  n2 <- a2$colnames
  nl1 <- a1$nonlinear
  nl2 <- a2$nonlinear
  as1 <- a1$assume.code
  as2 <- a2$assume.code

  l1 <- a1$len
  l2 <- a2$len
  if(any(nl1) & any(nl2)) nc <- l1+l2-1
  else nc <- l1*l2

  nr <- if(is.matrix(x1)) nrow(x1) else length(x1)

  x <- matrix(double(1),nrow=nr,ncol=nc)
  name <- character(nc)
  parms <- matrix(integer(1),nrow=2,ncol=nc+1)
  nonlinear <- logical(nc)

  k <- 0
  if(!is.factor(x1)) x1 <- as.matrix(x1)
  if(!is.factor(x2)) x2 <- as.matrix(x2)

  for(i in 1:l1) {
    if(as1==5 | as1==8) x1i <- unclass(x1)==(i+1)
    else x1i <- x1[,i]
    
    for(j in 1:l2) {
      ## Remove doubly nonlinear terms
      if(nl1[i] & nl2[j]) break
      
      k <- k + 1
      if(as2==5 | as2==8) x2j <- unclass(x2)==(j+1)
      else x2j <- x2[,j]
      
      x[,k] <- x1i * x2j
      name[k] <- paste(n1[i],"*",n2[j])
      parms[,k+1] <- c(nl1[i],nl2[j])
      nonlinear[k] <- nl1[i] | nl2[j]
    }
  }
  
  dimnames(x) <- list(NULL, name)
  attr(x,"ia") <- c(a1$name, a2$name)
  attr(x,"parms") <- parms
  attr(x,"nonlinear") <- nonlinear
  attr(x,"assume.code") <- 9
  attr(x,"name") <- paste(a1$name,"*",a2$name)
  attr(x,"label") <- attr(x,"name")
  attr(x,"colnames") <- name
  attr(x,"class") <- "rms"
  x
}


