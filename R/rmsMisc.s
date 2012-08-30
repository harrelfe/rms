#Miscellaneous functions to retrieve characteristics of design

#Function to get the number of intercepts in front of the slope coefficients
#ols - one intercept
#lrm - one or more intercepts (>1 if ordinal model)
#cph - no intercepts,   etc.

num.intercepts <- function(fit)
{
   nrp <- fit$non.slopes
   ## changed is.null(nrp) to below, fit$coefficients to fit$coef 14Aug01
   if(!length(nrp))
   {
	nm1 <- names(fit$coef)[1]  # 14Sep00
	nrp <- 1*(nm1=="Intercept" | nm1=="(Intercept)")
   }
   nrp
}

DesignAssign <- function(atr, non.slopes, Terms) {
  ## Given Design attributes and number of intercepts creates S-Plus
  ## format assign list (needed for R, intercept correction needed for
  ## S-Plus anyway).  If formula is given, names assign using
  ## terms(formul) term.labels, otherwise uses Design predictor names
  ## 23feb03: No, term.labels not useful if "." in formula
  ## formula argument no longer used

  ll <- if(missing(Terms)) atr$name else attr(Terms,'term.labels')
  if(!length(ll)) return(list())
  nv <- length(ll)
  params <- sapply(atr$nonlinear, length)  ## d.f. per predictor
  asc <- atr$assume.code
  assign <- list()
  j <- non.slopes + 1
  if(length(params)) for(i in 1:length(ll)) {
    if(asc[i]==8) next
    assign[[ll[i]]] <- j:(j+params[i]-1)
    j <- j + params[i]
  }
  assign
}
  
#Function to return variance-covariance matrix, optionally deleting
#rows and columns corresponding to parameters such as scale parameters
#in parametric survival models

vcov.lrm <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)  # for fastbw etc.
vcov.ols <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)
vcov.cph <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)
vcov.psm <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)

vcov.rms <- function(object, regcoef.only=TRUE, ...)
  {
    cov <- object$var
    if(!length(cov)) stop("fit does not have variance-covariance matrix")
    if(regcoef.only)
      {
        p <- length(object$coefficients)
        cov <- cov[1:p, 1:p, drop=FALSE]
      }
    cov
  }



## Functions for Out Of Sample computation of -2 log likelihood
## evaluated at parameter estimates of a given fit

oos.loglik <- function(fit, ...) UseMethod("oos.loglik")

oos.loglik.ols <- function(fit, lp, y, ...) {
  sigma2 <- sum(fit$residuals^2)/length(fit$residuals)
  if(missing(lp)) {
	n <- length(fit$residuals)
	n*logb(2*pi*sigma2)+n
  } else {
	s <- !is.na(lp + y)
	lp <- lp[s]; y <- y[s]
	n <- length(lp)
	sse <- sum((y - lp)^2)
	n*logb(2*pi*sigma2) + sse/sigma2
  }
}

oos.loglik.lrm <- function(fit, lp, y, ...) {
  if(missing(lp)) return(fit$deviance[length(fit$deviance)])
  ns <- fit$non.slopes
  if(ns > 1) stop('ordinal y case not implemented')
  y <- as.integer(as.category(y)) - 1
  s <- !is.na(lp + y)
  lp <- lp[s];  y <- y[s]
  p <- plogis(lp)
  -2*sum(ifelse(y==1, logb(p), logb(1-p)))
}
  
oos.loglik.cph <- function(fit, lp, y, ...) {
  if(missing(lp)) return(-2*fit$loglik[2])
  else stop('not implemented for cph models')
}

oos.loglik.psm <- function(fit, lp, y, ...) {
  if(missing(lp)) return(-2*fit$loglik[2])
  else stop('not implemented for psm models')
}

oos.loglik.Glm <- function(fit, lp, y, ...)
  if(missing(lp)) deviance(fit) else
  glm.fit(x=NULL, y=as.vector(y), offset=lp, family=fit$family)$deviance

  
#Function to retrieve limits and values, from fit (if they are there)
#or from a datadist object.  If need.all=F and input is coming from datadist,
#insert columns with NAs for variables not defined
#at is attr(fit$terms,"Design") (now fit$Design)

Getlim <- function(at, allow.null=FALSE, need.all=TRUE)
{
nam <- at$name[at$assume!="interaction"]
limits <- at$limits
values <- at$values

XDATADIST <- .Options$datadist
X <- lims <- vals <- NULL
if(!is.null(XDATADIST) && exists(XDATADIST))
  {
    X <- eval(as.name(XDATADIST))
    lims <- X$limits
    if(is.null(lims)) stop(paste("options(datadist=",XDATADIST,
                                 ") not created with datadist"))
    vals <- X$values
  }

if((length(X)+length(limits))==0) {
  if(allow.null) {
    lims <- list()
    for(nn in nam) lims[[nn]] <- rep(NA,7)
    lims <- structure(lims, class="data.frame", 
      row.names=c("Low:effect","Adjust to", "High:effect", "Low:prediction",
		  "High:prediction","Low","High"))
    return(list(limits=lims, values=values))
  }
  stop("no datadist in effect now or during model fit")
}

na <- if(length(limits))
  sapply(limits, function(x) all(is.na(x))) else rep(TRUE, length(nam))
if(length(lims) && any(na)) for(n in nam[na]) { #if() assumes NA stored in fit
						# for missing vars
  z <- limits[[n]]
  u <- if(match(n, names(lims), 0) > 0) lims[[n]] else NULL
  # This requires exact name match, not substring match
  if(is.null(u)) {
    if(need.all) stop(paste("variable",n,
	"does not have limits defined in fit or with datadist"))
    else limits[[n]] <- rep(NA,7)    # Added 28 Jul 94
  }
  else limits[[n]] <- u
}
limits <- structure(limits, class="data.frame", 
   row.names=c("Low:effect","Adjust to", "High:effect", "Low:prediction",
		"High:prediction","Low","High"))

if(length(vals)) values <- c(values, 
	vals[match(names(vals),nam,0)>0 & match(names(vals),names(values),0)==0]
	)   # add in values from datadist corresponding to vars in model
            # not already defined for model

list(limits=limits, values=values)
}

#Function to return limits for an individual variable, given an object
#created by Getlim

Getlimi <- function(name, Limval, need.all=TRUE)
{
   lim <- if(match(name, names(Limval$limits), 0) > 0) 
     Limval$limits[[name]] else NULL
   if(is.null(Limval) || is.null(lim) || all(is.na(lim))) {
      if(need.all) stop(paste("no limits defined by datadist for variable",
			name))
      return(rep(NA,7))
   }
lim
}

#Function to return a list whose ith element contains indexes
#of all predictors related, indirectly or directly, to predictor i
#Predictor i and j are related indirectly if they are related to
#any predictors that interact
#Set type="direct" to only include factors interacting with i
#This function is used by nomogram.

related.predictors <- function(at, type=c("all","direct"))
{
  type <- match.arg(type)
  f <- sum(at$assume.code < 9)
  if(any(at$assume.code == 10)) stop("does not work with matrix factors")
  ia <- at$interactions
  x <- rep(NA,f)
  names(x) <- at$name[at$assume.code < 9]
  mode(x) <- "list"
  if(length(ia)==0)
    {
      for(i in 1:f) x[[i]] <- integer(0)
      return(x)
    }
  for(i in 1:f)
    {
      r <- integer(0)
      for(j in 1:ncol(ia))
        {
          w <- ia[,j]
          if(any(w==i)) r <- c(r, w[w>0 & w!=i])
        }
      x[[i]] <- r
    }
  if(type=="direct") return(x)
  
  while(TRUE)
    {
      bigger <- FALSE
      for(j in 1:f)
        {
          xj <- x[[j]]
          y <- unlist(x[xj])
          y <- y[y != j]
          new <- unique(c(y, xj))
          bigger <- bigger | length(new) > length(xj)
          x[[j]] <- new
        }
      if(!bigger) break
    }
  x
}

#Function like related.predictors(..., type='all') but with new
# "super" predictors created by combining all indirected related
# (through interactions) predictors into a vector of predictor numbers
# with a new name formed from combining all related original names

combineRelatedPredictors <- function(at)
  {
    nam <- at$name
    r <- related.predictors(at)
    newnames <- newnamesia <- components <- list()
    pused <- rep(FALSE, length(nam))
    k <- 0
    for(i in (1:length(nam))[at$assume.code != 9])
      {
        if(!pused[i])
          {
            comp <- i
            nn   <- nam[i]
            ri   <- r[[i]]
            ianames <- character(0)
            ic <- interactions.containing(at, i)
            if(length(ic))
              {
                comp <- c(comp, ic)
                ianames <- nam[ic]
              }
            if(length(ri))
              {
                comp <- c(comp, ri)
                nn   <- c(nn,   nam[ri])
                for(j in ri)
                  {
                    pused[j] <- TRUE
                    ic <- interactions.containing(at, j)
                    if(length(ic))
                      {
                        comp <- c(comp, ic)
                        ianames <- c(ianames, nam[ic])
                      }
                  }
              }
            k <- k + 1
            components[[k]] <- unique(comp)
            newnames[[k]]   <- unique(nn)
            newnamesia[[k]] <- unique(c(nn, ianames))
          }
      }
    list(names=newnames, namesia=newnamesia, components=components)
  }
    

#Function to list all interaction term numbers that include predictor
#pred as one of the interaction components

interactions.containing <- function(at, pred) {
ia <- at$interactions
if(length(ia)==0) return(NULL)
name <- at$name
parms <- at$parms
ic <- NULL
for(i in (1:length(at$assume.code))[at$assume.code==9]) {
    terms.involved <- parms[[name[i]]][,1]
    if(any(terms.involved==pred)) ic <- c(ic, i)
}
ic
}

#Function to return a vector of logical values corresponding to
#non-intercepts, indicating if the parameter is one of the following types:
# term.order  Meaning
# ----------  -----------------
#     1       all parameters
#     2       all nonlinear or interaction parameters
#     3       all nonlinear parameters (main effects or interactions)
#     4       all interaction parameters
#     5       all nonlinear interaction parameters

param.order <- function(at, term.order) {	#at=Design attributes
if(term.order==1) return(rep(TRUE,length(at$colnames)))
nonlin <- unlist(at$nonlinear[at$name[at$assume!="strata"]]) # omit strat
ia <- NULL
for(i in (1:length(at$name))[at$assume!="strata"])
  ia <- c(ia, rep(at$assume[i]=="interaction",length(at$nonlinear[[i]])))
if(term.order==5) nonlin & ia else if(term.order==4) ia else
if(term.order==3) nonlin else nonlin | ia
}


#	rms.levels
#		Make each variable in an input data frame that is a
#		factor variable in the model be a factor variable with
#		the levels that were used in the model.  This is primarily
#		so that row insertion will work right with <-[.data.frame
#	
#at=Design attributes

rms.levels <- function(df, at)
{
  ac <- at$assume.code
  for(nn in names(df))
    {
      j <- match(nn, at$name, 0)
      if(j>0)
        {
          if((ac[j]==5 | ac[j]==8) & length(lev <- at$parms[[nn]]))
            df[[nn]] <- factor(df[[nn]], lev)
        }
    }
  df
}


#Function to return a default penalty matrix for penalized MLE,
#according to the design attributes and a design matrix X

Penalty.matrix <- function(at, X)
{
  d1 <- dimnames(X)[[2]][1]
  if(d1=='Intercept' || d1=='(Intercept)') X <- X[,-1,drop=FALSE]
  
  d <- dim(X)
  n <- d[1]; p <- d[2]
  center <- as.vector(rep(1/n,n) %*% X)   # see scale() function
  v <- as.vector(rep(1/(n-1),n) %*%
                 (X - rep(center,rep(n,p)))^2)
  
  pen <- if(p==1) as.matrix(v) else as.matrix(diag(v))    
  ## works even if X one column

  is <- 1
  ac <- at$assume
  for(i in (1:length(at$name))[ac!="strata"])
    {
      len <- length(at$nonlinear[[i]])
      ie <- is + len - 1
      if(ac[i] == "category") pen[is:ie,is:ie] <- diag(len) - 1/(len+1)
      is <- ie+1
    }
  pen
}

#Function to take as input a penalty specification of the form
#penalty=constant or penalty=list(simple=,nonlinear=,interaction=,
#nonlinear.interaction=) where higher order terms in the latter notation
#may be omitted, in which case their penalty factors are taken from lower-
#ordered terms.  Returns a new penalty object in full list form along
#with a full vector of penalty factors corresponding to the elements
#in regression coefficient vectors to be estimated

Penalty.setup <- function(at, penalty)
{
  if(!is.list(penalty))
    penalty <- list(simple=penalty, nonlinear=penalty,
                    interaction=penalty, nonlinear.interaction=penalty)
  tsimple <- penalty$simple
  if(!length(tsimple)) tsimple <- 0
  tnonlinear <- penalty$nonlinear
  if(!length(tnonlinear)) tnonlinear <- tsimple
  tinteraction <- penalty$interaction
  if(!length(tinteraction)) tinteraction <- tnonlinear
  tnonlinear.interaction <- penalty$nonlinear.interaction
  if(!length(tnonlinear.interaction)) tnonlinear.interaction <- tinteraction
  
  nonlin <- unlist(at$nonlinear[at$name[at$assume!='strata']])
  ia <- NULL
  for(i in (1:length(at$name))[at$assume!='strata'])
    ia <- c(ia, rep(at$assume[i]=='interaction',length(at$nonlinear[[i]])))
  nonlin.ia <- nonlin & ia
  nonlin[nonlin.ia] <- FALSE
  ia[nonlin.ia] <- FALSE
  simple <- rep(TRUE, length(ia))
  simple[nonlin | ia | nonlin.ia] <- FALSE
  penfact <- tsimple*simple + tnonlinear*nonlin + tinteraction*ia +
    tnonlinear.interaction*nonlin.ia
  list(penalty=list(simple=tsimple, nonlinear=tnonlinear,
         interaction=tinteraction,nonlinear.interaction=tnonlinear.interaction),
       multiplier=penfact)
}

#Function to do likelihood ratio tests from two models that are
# (1) nested and (2) have 'Model L.R.' components of the stats
# component of the fit objects
# For models with scale parameters, it is also assumed that the
# scale estimate for the sub-model was fixed at that from the larger model

lrtest <- function(fit1, fit2)
{
  if(length(fit1$fail) && fit1$fail)
    stop('fit1 had failed')
  if(length(fit2$fail) && fit2$fail)
    stop('fit2 had failed')
  
  s1 <- fit1$stats
  s2 <- fit2$stats
  
  if(!length(s1))
    s1 <- c('Model L.R.'=fit1$null.deviance - fit1$deviance,
            'd.f.'=fit1$rank - (any(names(coef(fit1))=='(Intercept)')))
  if(!length(s2))
    s2 <- c('Model L.R.'=fit2$null.deviance - fit2$deviance,
            'd.f.'=fit2$rank - (any(names(coef(fit2))=='(Intercept)')))
  
  chisq1 <- s1['Model L.R.']
  chisq2 <- s2['Model L.R.']
  if(length(chisq1)==0 || length(chisq2)==2) 
    stop('fits do not have stats component with "Model L.R." or deviance component')
  df1 <- s1['d.f.']
  df2 <- s2['d.f.']
  if(df1==df2) stop('models are not nested')

  lp1 <- length(fit1$parms);  lp2 <- length(fit2$parms)
  if(lp1 != lp2) warning('fits do not have same number of scale parameters') else 
  if(lp1 == 1 && abs(fit1$parms-fit2$parms)>1e-6)
    warning('fits do not have same values of scale parameters.\nConsider fixing the scale parameter for the reduced model to that from the larger model.')

  chisq <- abs(chisq1-chisq2)
  dof   <- abs(df1-df2)
  p     <- 1-pchisq(chisq,dof)

  r     <- c(chisq,dof,p)
  names(r) <- c('L.R. Chisq','d.f.','P')
  structure(list(stats=r,
                 formula1=formula(fit1),
                 formula2=formula(fit2)),
            class='lrtest')
}

print.lrtest <- function(x, ...)
{
  f1 <- x$formula1
  f2 <- x$formula2
  attributes(f1) <- NULL
  attributes(f2) <- NULL
  cat('\nModel 1: '); print(f1)
  cat('Model 2: '); print(f2); cat('\n')
  print(x$stats)
  cat('\n')
  invisible()
}


Newlabels <- function(fit, ...) UseMethod('Newlabels')

Newlabels.rms <- function(fit, labels, ...)
{
  at <- fit$Design
  nam <- names(labels)
  if(length(nam)==0)
    {
      if(length(labels)!=length(at$name))
        stop('labels is not a named vector and its length is not equal to the number of variables in the fit')
      nam <- at$name
    } 
  i <- match(nam, at$name, nomatch=0)
  
  if(any(i==0))
    {
      warning(paste('the following variables were not in the fit and are ignored:\n',
                    paste(nam[i==0],collapse=' ')))
      labels <- labels[i>0]
      i <- i[i>0]
    }
  
  at$label[i] <- labels
  
  fit$Design <- at
  fit
}

Newlevels <- function(fit, ...) UseMethod('Newlevels')

Newlevels.rms <- function(fit, levels, ...)
{
  at <- fit$Design
  nam <- names(levels)
  if(length(nam)==0) stop('levels must have names')

  i <- match(nam, at$name, nomatch=0)
  
  if(any(i==0))
    {
      warning(paste('the following variables were not in the fit and are ignored:\n',
                    paste(nam[i==0],collapse=' ')))
      nam <- nam[i>0]
    }
  
  for(n in nam)
    {
      prm <- at$parms[[n]]
      if(length(prm)!=length(levels[[n]]))
        stop(paste('new levels for variable',
                   n,'has the wrong length'))

      levs <- levels[[n]]
      if(length(at$values[[n]])) at$values[[n]] <- levs
      if(length(at$limits))
        {
          m <- match(at$limits[[n]], at$parms[[n]])
          if(is.category(at$limits[[n]]))
            attr(at$limits[[n]],'levels') <- levs
          else
            at$limits[[n]] <- levs[m]
        }
      at$parms[[n]] <- levs
    }
  
  fit$Design <- at
  fit
}

rmsFit <- function(fit)
{
  cl <- oldClass(fit)
  if(cl[1]=='rms') return(fit)
  fit$fitFunction <- cl
  oldClass(fit) <- 'rms'
  fit
}

print.rms <- function(x, ...)
{
  fitter <- x$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(x) <- fitter[1]
  print(x, ...)
}

residuals.rms <- function(object, ...)
{
  fitter <- object$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(object) <- fitter[1]
  residuals(object, ...)
}

validate.rms <- function(fit, ...)
{
  fitter <- fit$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(fit) <- fitter[1]
  validate(fit, ...)
}

calibrate.rms <- function(fit, ...)
{
  fitter <- fit$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(fit) <- fitter[1]
  calibrate(fit, ...)
}

Survival.rms <- function(object, ...)
{
  fitter <- object$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(object) <- fitter[1]
  Survival(object, ...)
}

Quantile.rms <- function(object, ...)
{
  fitter <- object$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(object) <- fitter[1]
  Quantile(object, ...)
}

Mean.rms <- function(object, ...)
{
  fitter <- object$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(object) <- fitter[1]
  Mean(object, ...)
}

Hazard.rms <- function(object, ...)
{
  fitter <- object$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(object) <- fitter[1]
  Hazard(object, ...)
}

latex.rms <-
  function(object, title,
           file=paste(first.word(deparse(substitute(object))),
             'tex',sep='.'), ...)
{
  fitter <- object$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  oldClass(object) <- fitter[1]
  ## Need to brute-force dispatch because of SV4 problem in latex in Hmisc
  if(existsFunction(p <- paste('latex',fitter[1],sep='.')))
    do.call(p, list(object, file=file, ...))
  else
    latexrms(object, file=file, ...)
}

survest.rms <- function(fit, ...)
{
  fitter <- fit$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  f <- paste('survest',fitter[1],sep='.')
  do.call(f, list(fit,...))
}


oos.loglik.rms <- function(fit, ...)
{
  fitter <- fit$fitFunction
  if(!length(fitter))
    stop("fit's main class is 'rms' but no fitFunction element is present")
  f <- paste('oos.loglik',fitter[1],sep='.')
  do.call(f, list(fit,...))
}

#getOldDesign <- function(fit) {
#  at <- attr(fit$terms,'Design')
#  if(is.null(at))
#    stop('fit was not created by a Design library fitting function')
#  at
#}

univarLR <- function(fit)
{
  ## Computes all univariable LR chi-square statistics
  w <- as.character(attr(fit$terms,'variables'))
  w <- w[-1]
  p <- length(w)-1
  stats <- P <- double(p)
  dof <- nobs <- integer(p)
  for(i in 1:p)
    {
      stat <- update(fit, as.formula(paste(w[1],w[i+1],sep='~')))$stats
      stats[i] <- stat['Model L.R.']
      dof[i]   <- stat['d.f.']
      P[i]     <- stat['P']
      nobs[i]  <- stat['Obs']
    }
  data.frame(LR=stats, 'd.f.'=dof, P=P, N=nobs,
             row.names=w[-1], check.names=FALSE)
}

vif <- function(fit)
{
  v <- vcov(fit, regcoef.only=TRUE)
  nam <- dimnames(v)[[1]]
  ns <- num.intercepts(fit)
  if(ns>0) {
    v <- v[-(1:ns),-(1:ns),drop=FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

## Returns a list such that variables with no = after them get the value NA
## For handling ... arguments to Predict, summary, nomogram, gendata,
## survplot.rms, ...
rmsArgs <- function(.object, envir=parent.frame(2))
  {
    if(length(.object) < 2) return(NULL)
    .names <- names(.object)[-1]
    ## See if no variables given with = after their names
    if(!length(.names)) .names <- rep('', length(.object)-1)
    .n <- length(.names)
    .vars  <- sapply(.object, as.character)[-1]
    .res <- vector('list', .n)
    for(.i in 1:.n)
      {
        if(.names[.i] == '')
          {
            .names[.i] <- .vars[.i]
            .res[[.i]] <- NA
          }
        else .res[[.i]] <- eval(.object[[.i+1]], envir=envir)
      }
    names(.res) <- .names
    .res
  }

## General function to print model fit objects using latex or regular
## print (the default)

prModFit <- function(x, title, w, digits=4, coefs=TRUE,
                     latex=FALSE, lines.page=40, long=TRUE, needspace, ...)
  {
    bverb <- function() if(latex) cat('\\begin{verbatim}\n')
    everb <- function() if(latex) cat('\\end{verbatim}\n')
    skipt  <- function(n=1, latex=FALSE)
      {
        if(n==0) return()
        if(latex) paste('\n\\vspace{', n, 'ex}\n\n', sep='')
        else paste(rep('\n', n), collapse='')
      }
    catl  <- function(x, skip=1, bold=FALSE, verb=FALSE, pre=0,
                      center=TRUE, indent=FALSE)
      if(latex)
        {
          if(verb)
            cat(paste('\\begin{verbatim}\n', skipt(pre), sep=''),
                x,
                skipt(skip),
                '\\end{verbatim}\n',
                sep='')
          else
            cat(skipt(pre, latex=TRUE),
                if(center) '\\centerline{'
                else if(!indent) '\\noindent ',
                if(bold) '\\textbf{',
                x,
                if(bold) '}',
                if(center) '}',
                skipt(skip, latex=TRUE),
                sep='')
        }
      else
        cat(skipt(pre), x, '\n', skipt(skip), sep='')
    
    
    latexVector <- function(x, ...)
      cat(latexTabular(t(x), helvetica=FALSE, ...),'\n', sep='')
    
    if(length(x$fail) && x$fail)
      {
        catl('Model Did Not Converge.  No summary provided.', bold=TRUE, pre=1)
        return()
      }

    if(!missing(needspace) && latex)
      cat('\\Needspace{', needspace, '}\n', sep='')
    if(title != '') catl(title, pre=1, bold=TRUE)

    if(long)
      {
        bverb()
        dput(x$call)
        cat('\n')
        everb()
      }

    for(z in w)
      {
        type <- z$type
        obj  <- z[[2]]
        titl <- z$title
        tex  <- z$tex
        if(!length(tex)) tex <- FALSE
        if(type == 'naprint.delete' && latex) {
          type <- 'latex.naprint.delete'
          tex <- TRUE
        }

        preskip <- z$preskip
        if(!length(preskip)) preskip <- 0
        if(!tex && length(titl)) catl(titl, pre=preskip, skip=1)
        if(type == 'stats')
          {
            prStats(obj[[1]], obj[[2]], latex=latex)
            if(!latex) cat('\n')
          }
        else if(type == 'coefmatrix')
          {
            if(coefs)
              {
                errordf <- obj$errordf
                beta <- obj$coef
                se   <- obj$se
                Z    <- beta/se
                P    <- if(length(errordf)) 2*(1 - pt(abs(Z), errordf)) else
                        1 - pchisq(Z^2, 1)
                pad <- function(x)
                  if(latex) paste('~', x, '~', sep='') else x
                U    <- cbind('\\textrm{~Coef~}' =
                              pad(formatNP(beta, digits, latex=latex)),
                              '\\textrm{~S.E.~}' =
                              pad(formatNP(se,   digits, latex=latex)),
                              '\\textrm{Wald~} Z'  =
                              formatNP(Z,    2, latex=latex),
                              '\\textrm{Pr}(>|Z|)' =
                              formatNP(P, 4, latex=latex, pvalue=TRUE))
                if(!latex)
                  colnames(U) <- c('Coef', 'S.E.', 'Wald Z', 'Pr(>|Z|)')
                if(length(errordf))
                  colnames(U)[3:4] <- if(latex) c('t', '\\textrm{Pr}(>|t|)') else
                                                c('t',   'Pr(>|t|)')
                rownames(U) <- names(beta)
                if(length(obj$aux))
                  {
                    U <- cbind(U, formatNP(obj$aux, digits, latex=latex))
                    colnames(U)[ncol(U)] <- obj$auxname
                  }
                if(latex)
                  {
                    cat(skipt(1, latex=TRUE))
                    rownames(U) <- latexTranslate(names(beta))
                    if(is.numeric(coefs))
                      {
                        U <- U[1:coefs,,drop=FALSE]
                        U <- rbind(U, rep('', ncol(U)))
                        rownames(U)[nrow(U)] <- '\\dots'
                      }
                    if(!missing(needspace) && latex)
                      cat('\\Needspace{', needspace, '}\n', sep='')
                    latex(U, file='', first.hline.double=FALSE,
                          table=FALSE, longtable=TRUE,
                          lines.page=lines.page,
                          col.just=rep('r',4), rowlabel='',
                          math.col.names=TRUE)
                  }
                else
                  {
                    if(is.numeric(coefs))
                      {
                        U <- U[1:coefs,,drop=FALSE]
                        U <- rbind(U, rep('', ncol(U)))
                        rownames(U)[nrow(U)] <- '. . .'
                      }
                    print(U, quote=FALSE)
                    cat('\n')
                  }
              }
          }
        else {
          if(tex)
            {
              cat('\\begin{center}\n')
              if(length(titl)) cat(titl, '\n\n')
            }
          else
            {
              bverb()
              cat(skipt(preskip, latex=tex))
            }
          do.call(type, obj)
          ## unlike do.call, eval(call(...)) dispatches on class of ...
          if(tex) cat('\\end{center}\n')
          else
            {
              cat('\n')
              everb()
            }
        }
      }
    cat('\n')
  }

latex.naprint.delete <- function(object, ...) {
  lg <- length(g <- object$nmiss)
  if(lg) {
    cat("Frequencies of Missing Values Due to Each Variable\n\n\\smallskip\n\n")
    if(sum(g > 0) < 4) {
      cat('\\begin{verbatim}\n')
      print(g)
      cat('\\end{verbatim}\n')
    } else {
      z <- latexDotchart(g, names(g), auxdata=g, auxtitle='N',
                         w=3, h=min(max(2.5*lg/20, 1), 8))
      cat(z, sep='\n')
    }
    cat("\n")
  }
  
  if(length(g <- object$na.detail.response)) {
    cat("\nStatistics on Response by Missing/Non-Missing Status of Predictors\n\n")
    print(oldUnclass(g))
    cat("\n")           
  }
  
  invisible()
}
                         
    
## Function to print model fit statistics
## Example:
#prStats(list('Observations', c('Log','Likelihood'),
#            c('Rank','Measures'),
#            c('Mean |difference|','Measures')),
#       list(c(N0=52, N1=48), c('max |deriv|'=1e-9,'-2 LL'=1332.23,c(NA,2)),
#            c(tau=-.75, Dxy=-.64, C=.743, 2),
#            c(g=1.25, gr=11.3, 2)))
## Note that when there is an unnamed element of w, it is assumed to be
## the number of digits to the right of the decimal place (recycling of
## elements is done if fewer elements are in this vector), causing
## round(, # digits) and format(..., nsmall=# digits).  Use NA to use
## format without nsmall and without rounding (useful for integers and for
## scientific notation)

prStats <- function(labels, w, latex=FALSE)
  {
    spaces <- function(n)
      if(n <= 0.5) '' else
    substring('                                                         ',
              1, floor(n))
    
    ## Find maximum width used for each column
    p <- length(labels)
    width <- numeric(p)
    for(i in 1:p)
      {
        width[i] <- max(nchar(labels[[i]]))
        u <- w[[i]]
        dig <- NA
        if(any(names(u)==''))
          {
            dig <- u[names(u)=='']
            u   <- u[names(u)!='']
          }
        lu <- length(u)
        dig <- rep(dig, length=lu)
        fu <- character(lu)
        for(j in 1:length(u))
          {
            dg <- dig[j]
            fu[j] <- if(is.na(dg)) format(u[j]) else
            if(dg < 0) formatNP(u[j], -dg, pvalue=TRUE, latex=latex) else
            formatNP(u[j], dg, latex=latex)
          }
        names(fu) <- names(u)
        w[[i]] <- fu
        for(j in 1:length(u))
          width[i] <- max(width[i],
                          1 + nchar(names(u))[j] + nchar(fu[j]))
      }
    if(latex)
      {
        cat('\\centerline{\\begin{tabular}{|', rep('c|',p), '}\\hline\n',
            sep='')
        if(sum(nchar(unlist(labels))) > 0)
          {
            maxl <- max(sapply(labels, length))
            for(i in 1:maxl)
              {
                lab <- sapply(labels, function(x) if(length(x) < i)'' else x[i])
                cat(paste(lab, collapse='&'), '\\\\ \n', sep='')
              }
            cat('\\hline\n')
          }
        maxl <- max(sapply(w, length))
        z <- matrix('', nrow=maxl, ncol=p)
        for(i in 1:p)
          {
            k <- latexTranslate(names(w[[i]]), greek=TRUE)
            k[k=='Dxy']   <- '$D_{xy}$'
            k[k=='LR chi2']  <- 'LR $\\chi^{2}$'
            k[k=='Score chi2'] <- 'Score $\\chi^{2}$'
            k[k=='Pr($>$ chi2)'] <- 'Pr$(>\\chi^{2})$'
            k[k=='$\\tau$-a'] <- '$\\tau_{a}$'
            k[k=='R2']    <- '$R^{2}$'
            k[k=='R2 adj'] <- '$R^{2}_{\\textrm{adj}}$'
            k[k=='C']     <- '$C$'
            k[k=='g']     <- '$g$'
            k[k=='gp']    <- '$g_{p}$'
            k[k=='gr']    <- '$g_{r}$'
            k[k=='max $|$deriv$|$'] <- 'max $|$deriv$|$~'
            z[1:length(k),i] <- paste(k, '~\\hfill ', w[[i]], sep='')
          }
        for(j in 1:maxl) cat(paste(z[j,], collapse='&'), '\\\\ \n', sep='')
        cat('\\hline\n')
        cat('\\end{tabular}}\n\n')
        return()
      }
    z <- labs <- character(0)
    for(i in 1:p)
      {
        wid <- width[i]
        lab <- labels[[i]]
        for(j in 1:length(lab))
          lab[j] <- paste(spaces((wid - nchar(lab[j]))/2), lab[j], sep='')
        labs <- c(labs, paste(lab, collapse='\n'))
        u   <- w[[i]]
        a <- ''
        for(i in 1:length(u))
          a <- paste(a, names(u)[i],
                     spaces(wid - nchar(u[i]) - nchar(names(u[i]))),
                     u[i],
                     if(i < length(u)) '\n', sep='')
        z <- c(z, a)
      }
    res <- rbind(labs, z)
    rownames(res) <- NULL
    print.char.matrix(res, vsep='', hsep='    ', csep='',
                      top.border=FALSE, left.border=FALSE)
  }

## reVector is used in conjunction with pstats
## Example:
# x <- c(a=1, b=2)
# c(A=x[1], B=x[2])
# reVector(A=x[1], B=x[2])
# reVector(A=x['a'], B=x['b'], C=x['c'])
reVector <- function(..., na.rm=TRUE)
  {
    d <- list(...)
    d <- d[sapply(d, function(x) !is.null(x))]
    x <- unlist(d)
    names(x) <- names(d)
    if(na.rm) x[!is.na(x)] else x
  }

formatNP <- function(x, digits=NULL, pvalue=FALSE, latex=FALSE)
  {
    f <- if(length(digits))
      format(round(x, digits), nsmall=digits, scientific=1) else
      format(x, scientific=1)
    sci <- grep('e', f)
    if(latex && length(sci)) f[sci] <- paste('$', latexSN(f[sci]), '$', sep='')
    if(!pvalue) return(f)
    if(!length(digits)) stop('must specify digits if pvalue=TRUE')
    s <- x < 10^(-digits)
    if(any(s))
      {
        w <- paste('0.', paste(rep('0', digits-1), collapse=''), '1', sep='')
        f[s] <- if(latex) paste('$<', w, '$', sep='') else
        paste('<', w, sep='')
      }
    f
  }

logLik.ols <- function(object, ...)
  stats:::logLik.lm(object)

logLik.rms <- function(object, ...)
  {
    dof <- unname(object$stats['d.f.'] + num.intercepts(object))
    if(inherits(object, 'psm')) dof <- dof + 1  # for sigma
    nobs <- nobs(object)
    w <- object$loglik
    if(length(w)) return(structure(w[length(w)], nobs=nobs, df=dof, class='logLik'))
    w <- object$deviance
    structure(-0.5*w[length(w)], nobs=nobs, df=dof, class='logLik')
  }

AIC.rms <- function(object, ..., k=2, type=c('loglik','chisq'))
  {
    type <- match.arg(type)
    if(type == 'loglik') return(AIC(logLik(object), k=k))
    stats <- object$stats
    dof <- stats['d.f.']
    unname(stats['Model L.R.'] - k * dof)
  }

nobs.rms <- function(object, ...)
  {
    st <- object$stats
    if(inherits(object,'Gls')) length(object$residuals)
    else if(any(names(st) == 'Obs')) unname(st['Obs'])
    else unname(st['n'])
  }
