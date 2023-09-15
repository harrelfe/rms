#Miscellaneous functions to retrieve characteristics of design


DesignAssign <- function(atr, non.slopes, Terms) {
  ## Given Design attributes and number of intercepts creates R
  ## format assign list.

  ll <- if(missing(Terms)) atr$name else attr(Terms,'term.labels')
  if(! length(ll)) return(list())
  nv     <- length(ll)
  params <- sapply(atr$nonlinear, length)  ## d.f. per predictor
  asc    <- atr$assume.code
  assign <- list()
  j <- non.slopes + 1
  if(length(params)) for(i in 1 : length(ll)) {
    if(asc[i] == 8) next
    assign[[ll[i]]] <- j : (j + params[i] - 1)
    j <- j + params[i]
  }
  assign
}
  
#Function to return variance-covariance matrix, optionally deleting
#rows and columns corresponding to parameters such as scale parameters
#in parametric survival models (if regcoef.only=TRUE)

vcov.lrm <- function(object, regcoef.only=TRUE, intercepts='all', ...) {
  if(length(intercepts) == 1 && is.character(intercepts) &&
     intercepts %nin% c('all', 'none'))
    stop('if character, intercepts must be "all" or "none"')

  if(!length(intercepts) ||
     (length(intercepts) == 1) && intercepts == 'all')
    return(vcov.rms(object, regcoef.only=regcoef.only, ...))
  ns <- num.intercepts(object)
  v <- object$var
  p <- ncol(v)
  nx <- p - ns
  if(intercepts == 'none') intercepts <- integer(0)
  i <- if(nx == 0) intercepts else c(intercepts, (ns+1):p)
  v[i, i, drop=FALSE]
}

vcov.ols <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)

vcov.cph <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)

vcov.psm <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)

vcov.orm <- function(object, regcoef.only=TRUE,
                     intercepts='mid', ...) {
  v <- object$var
  if(! length(intercepts)) return(v)
  li1 <- length(intercepts) == 1
  iat <- attr(v, 'intercepts')  # handle fit.mult.impute (?), robcov
  # robcov re-writes var object and uses all intercepts
  iref <- object$interceptRef
  if(is.numeric(intercepts) && li1 &&
     intercepts == iref) intercepts <- 'mid'
  if(! length(iat)) {
    if(li1 && intercepts == 'mid') {
      i <- c(iref, (num.intercepts(object, 'var') + 1) : nrow(v))
      return(object$var[i, i, drop=FALSE])
    }
    return(vcov.lrm(object, regcoef.only=regcoef.only,
                    intercepts=intercepts, ...))
  }
  
  if(li1 && intercepts == 'none')
    return(object$var[-(1 : length(iat)),
                      -(1 : length(iat)), drop=FALSE])
    if(li1 && intercepts == 'mid' && length(iat) == 1) return(object$var)
  
  iref <- object$interceptRef
  info <- object$info.matrix
  isbootcov <- length(object$boot.coef)
  ns <- num.intercepts(object)
  p  <- ncol(info)
  ns <- num.intercepts(object)
  nx <- p - ns
  scale <- attr(info, 'scale')
  name <- names(coef(object))
  if(length(scale) && (! is.character(intercepts) ||
                       (li1 && intercepts == 'all'))) {
    xbar  <- scale$mean
    xsd   <- scale$sd
    trans <- 
      rbind(cbind(diag(ns), matrix(0, nrow=ns, ncol=nx)),
            cbind(-matrix(rep(xbar / xsd, ns), ncol=ns),
                  diag(1 / as.vector(xsd))))
  }
  
                       
  if(li1 && is.character(intercepts)) {
    if(intercepts != 'mid' && isbootcov)
      stop('intercepts must be "mid" if object produced by bootcov')
      switch(intercepts,
             mid = return(object$var),
             all = {
               if(! length(scale)) {
                 v <- as.matrix(solve(info))
                 dimnames(v) <- list(name, name)
                 return(v)
               }
               kint <- num.intercepts(object)
               v <- t(trans) %*% as.matrix(solve(info)) %*% trans
               dimnames(v) <- list(name, name)
               return(v)
             },
             none= return(object$var[-1, -1, drop=FALSE]) )
  }
  if(isbootcov)
    stop('intercepts must be "mid" if object produced by bootcov')
  
  i <- if(nx == 0) intercepts else c(intercepts, (ns+1):p)
  v <- if(length(scale))
    (t(trans) %*% as.matrix(solve(info)) %*% trans)[i,i]
  else as.matrix(solve(info)[i,i])
  dimnames(v) <- list(name[i], name[i])
  v
}

vcov.rms <- function(object, regcoef.only=TRUE, intercepts='all', ...)
  {
    cov <- object$var
    if(!length(cov)) stop("fit does not have variance-covariance matrix")
    if(regcoef.only) {
      p <- length(object$coefficients)
      cov <- cov[1:p, 1:p, drop=FALSE]
    }
    if(length(intercepts) && intercepts == 'none') {
      ns <- num.intercepts(object)
      if(ns > 0) cov <- cov[-(1:ns), -(1:ns), drop=FALSE]
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
  y <- as.integer(as.factor(y)) - 1
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

Getlim <- function(at, allow.null=FALSE, need.all=TRUE) {
nam    <- at$name[at$assume!="interaction"]
limits <- at$limits
values <- at$values

XDATADIST <- .Options$datadist
X <- lims <- vals <- NULL
if(! is.null(XDATADIST)) {
  X <- if(inherits(XDATADIST, 'datadist')) XDATADIST
       else
         if(exists(XDATADIST)) eval(as.name(XDATADIST))
  if(! is.null(X)) {
    lims <- X$limits
    if(is.null(lims)) stop(paste("options(datadist=",XDATADIST,
                                 ") not created with datadist"))
    vals <- X$values
    }
  }

if((length(X) + length(limits)) == 0) {
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
  if(d1 %in% c('Intercept', '(Intercept)')) X <- X[, -1, drop=FALSE]
  
  d <- dim(X)
  n <- d[1]; p <- d[2]
  center <- as.vector(rep(1 / n, n) %*% X)   # see scale() function
  v <- as.vector(rep(1 / (n - 1), n) %*%
                 (X - rep(center, rep(n, p)))^2)
  
  pen <- if(p == 1) as.matrix(v) else as.matrix(diag(v))    
  ## works even if X one column

  is <- 1
  ac <- at$assume
  for(i in (1 : length(at$name))[ac != "strata"]) {
    len <- length(at$nonlinear[[i]])
    ie <- is + len - 1
    if(ac[i] == "category") pen[is : ie, is : ie] <- diag(len) - 1 / (len + 1)
    is <- ie + 1
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
          if(is.factor(at$limits[[n]]))
            attr(at$limits[[n]],'levels') <- levs
          else
            at$limits[[n]] <- levs[m]
        }
      at$parms[[n]] <- levs
    }
  
  fit$Design <- at
  fit
}

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

## General function to print model fit objects using latex, html, or regular
## print (the default)

prModFit <- function(x, title, w, digits=4, coefs=TRUE, footer=NULL,
                     lines.page=40, long=TRUE, needspace, subtitle=NULL, ...) {
  lang   <- prType()
  specs  <- markupSpecs[[lang]]
  transl <- switch(lang,
                   latex = latexTranslate,
                   html  = htmlTranslate,
                   plain = function(x) x)

#  cca  <- htmlSpecial('combiningcircumflexaccent')
  nbsp <- htmlSpecial('nbsp')
  gt   <- transl('>')
  vbar <- transl('|')
  chi2 <- specs$chisq()
  beta <- htmlGreek('beta')
  
  R <- character(0)

  bverb <- function() {
    switch(lang,
           html  = '<pre>',
           latex = '\\begin{verbatim}',
           plain = NULL)
    }
  
  everb <- function()
    switch(lang,
           html  = '</pre>',
           latex = '\\end{verbatim}',
           plain = NULL)
  
  skipt  <- function(n=1) {
    if(n==0) return(character(0))
    if(n == 1) return('')
    specs$lineskip(n)
  }

  catl  <- function(x, skip=1, bold=FALSE, verb=FALSE, pre=0,
                    center=FALSE, indent=FALSE) {
    if(lang == 'latex') {
      if(verb)
        c('\\begin{verbatim}', skipt(pre),
          x,
          skipt(skip),
          '\\end{verbatim}')
      else
        c(skipt(pre),
          paste0(
            if(center) '\\centerline{'
            else if(!indent) '\\noindent ',
            if(bold) '\\textbf{',
            x,
            if(bold) '}',
            if(center) '}'),
          skipt(skip))
    } else if(lang == 'html') {
      if(verb)
        c('<pre>', skipt(pre),
          x,
          skipt(skip),
          '</pre>')
      else
        c(skipt(pre),
          paste0(if(center) '<div align=center>' else '<p>',
                 if(bold) '<strong>',
                 x,
                 if(bold) '</strong>',
                 if(center) '</div>' else '</p>'),
          skipt(skip))
    }
    else c(paste0(skipt(pre), x), skipt(skip))
  }
  latexVector <- function(x, ...)
    latexTabular(t(x), helvetica=FALSE, ...)
  
  if(length(x$fail) && x$fail) {
    return(catl('Model Did Not Converge.  No summary provided.',
                bold=TRUE, pre=1, verb=TRUE))
  }

  R <- character(0)
  
  if(! missing(needspace) && lang == 'latex')
    R <- paste0('\\Needspace{', needspace, '}')

  lsub <- length(subtitle)
  if(title != '') R <- c(R, catl(title, pre=1, bold=TRUE,
                                 skip=1))
  ## was skip=if(lsub) 0 else 1
  if(lsub)
    for(i in lsub) R <- c(R, catl(subtitle[i], bold=FALSE, pre=1))
  
  if(long) {
    R <- c(R, bverb(), deparse(x$call), everb(), '')
    ## dput(x$call) didn't work with rmarkdown because dput has no append=
  }

  for(z in w) {
    type <- z$type
    obj  <- z[[2]]
    titl <- z$title
    tex  <- z$tex
    if(! length(tex)) tex <- FALSE
    if(type == 'naprint.delete') {
      if(lang == 'latex') {
        type <- 'latex.naprint.delete'
        tex <- TRUE
      }
      if(lang == 'html') type <- 'html.naprint.delete'
    }
    
    preskip <- z$preskip
    if(! length(preskip)) preskip <- 0
    if(! tex && length(titl)) R <- c(R, '', catl(titl, pre=preskip, skip=1))
    if(type == 'stats') {
      R <- c(R, prStats(obj[[1]], obj[[2]], lang=lang))
    } else if(type == 'coefmatrix') {
      if(coefs) {
        pad <- function(x)
          switch(lang, 
                 latex = paste0('~', x, '~'),
                 html  = paste0(nbsp, x),
                 plain  = x)
        betan <- switch(lang,
                        plain = 'Beta',
                        html  = htmlGreek('beta'),
                        latex = '$\\hat{\\beta}$')
                        
        B   <- obj$bayes
        if(length(B)) {
          U <- matrix('', nrow=nrow(B), ncol=ncol(B))
          for(i in 1:ncol(B)) {
            dig <- if(colnames(B)[i] == 'Symmetry') 2 else digits
            U[, i] <- pad(formatNP(B[, i], dig, lang=lang))
          }
          pn <- switch(lang, plain='Pr(Beta>0)',
                       html = paste0('Pr(', betan, transl('>'), '0)'),
                       latex = 'Pr$(\\beta>0)$')
          coltrans <- c(Mean     = paste('Mean', betan),
                        Median   = paste('Median', betan),
                        Mode     = paste('Mode', betan),
                        SE       = 'S.E.',
                        Lower    = 'Lower',
                        Upper    = 'Upper',
                        P        = pn,
                        Symmetry = 'Symmetry')
          colnames(U) <- coltrans[colnames(B)]
          rownames(U) <- rownames(B)
          betanames   <- rownames(B)
        }
        else  {
        errordf <- obj$errordf
        beta <- obj$coef
        betanames <- names(beta)
        se   <- obj$se
        Z    <- beta / se
        P    <- if(length(errordf)) 2 * (1 - pt(abs(Z), errordf))
                else
                  1 - pchisq(Z ^ 2, 1)

        U    <- cbind('Coef' =
                        pad(formatNP(beta, digits, lang=lang)),
                      'S.E.' =
                        pad(formatNP(se,   digits, lang=lang)),
                      'Wald Z'  =
                        formatNP(Z,    2, lang=lang),
                      'Pr(>|Z|)' =
                        formatNP(P, 4, lang=lang, pvalue=TRUE))
        if(lang == 'latex')
          colnames(U) <- c('$\\hat{\\beta}$', 'S.E.', 'Wald $Z$',
                           'Pr$(>|Z|)$')
        else
          if(lang == 'html')
            colnames(U) <- c(htmlGreek('beta'),   # did have cca
                             'S.E.', 'Wald <i>Z</i>',
                             paste0('Pr(', gt, vbar, '<i>Z</i>', vbar, ')'))
        if(length(errordf))
          colnames(U)[3:4] <-
            switch(lang,
                   latex = c('$t$', 'Pr$(>|t|)$'),
                   html  = c('<i>t</i>', paste0('Pr(', gt, vbar, '<i>t</i>',
                                                vbar, ')')),
                   plain = c('t',   'Pr(>|t|)') )

        rownames(U) <- betanames

        if(length(obj$aux)) {
          U <- cbind(U, formatNP(obj$aux, digits, lang=lang))
          colnames(U)[ncol(U)] <- obj$auxname
        }
        }
        if(lang %in% c('latex', 'html')) {
          R <- c(R, skipt(1))
          rownames(U) <- transl(betanames)
          
          if(is.numeric(coefs)) {
            U <- U[1:coefs,,drop=FALSE]
            U <- rbind(U, rep('', ncol(U)))
            rownames(U)[nrow(U)] <- if(lang == 'html') '&hellip;' else '\\dots'
          }
          ## Translate interaction symbol (*) to times symbol
          rownames(U) <- gsub('*', specs$times, rownames(U), fixed=TRUE)
 
          if(! missing(needspace) && lang == 'latex')
            R <- c(R, paste0('\\Needspace{', needspace, '}'))

          if(lang == 'latex') 
            R <- c(R,   # was capture.output(latex())
                   capture.output(latex(U, file='',
                                        first.hline.double=FALSE,
                                        table=FALSE, longtable=TRUE,
                                        lines.page=lines.page,
                                        col.just=rep('r',ncol(U)), rowlabel='',
                                        already.math.col.names=TRUE,
                                        append=TRUE)))
          else {
            al <- paste(rep('r', ncol(U)), collapse='')
            R <- c(R, as.character(
                        htmlTable::htmlTable(U,
                                             css.cell = 'min-width: 7em;',
                                             align=al, align.header=al,
                                             # rowlabel='',
                                             escape.html=FALSE)))
            }
        } else {
          if(is.numeric(coefs)) {
            U <- U[1:coefs,,drop=FALSE]
            U <- rbind(U, rep('', ncol(U)))
            rownames(U)[nrow(U)] <- '. . .'
          }
          R <- c(R, '', capture.output(print(U, quote=FALSE)), '')
        }
      }   ## end if(coefs)
    }     ## end coefmatrix
    else {
      if(tex) {    ### ??? how does this apply to html?
        R <- c(R, '\\begin{center}',
               if(length(titl)) c(titl, '\n'))
      } else {
        R <- c(R,  skipt(preskip))
      }
      R <- c(R,
             if(type == 'html.naprint.delete')
               do.call(type, obj)
             else
               if(type == 'latex.naprint.delete')
                 capture.output(do.call(type,
                                        c(obj, list(file=''))))
             else
               if(type == 'print')
                 c(bverb(),
                   capture.output(do.call(type,
                                          c(obj, list(quote=FALSE)))), everb())
             else
               do.call(type, obj),
             ## unlike do.call, eval(call(...)) dispatches on class of ...
             if(tex) '\\end{center}' else ''
      )
    }
  }
  if(length(footer))
    R <- c(R, paste(specs$smallskip, transl(footer)))

  if(getOption('rmsdebug', FALSE))
    cat(R, sep='\n', append=TRUE, file='/tmp/rmsdebug.txt')

  switch(lang,
         html  = rendHTML(R),
         latex = cat(R, sep='\n'),
         plain = cat(R, sep='\n'))
}

latex.naprint.delete <- function(object, file='', append=TRUE, ...) {
  lg <- length(g <- object$nmiss)
  if(file != '') sink(file, append=append)
  if(lg) {
    cat("Frequencies of Missing Values Due to Each Variable\n\n\\smallskip\n\n")
    if(sum(g > 0) < 4) {
      cat('\\begin{verbatim}\n')
      print(g)
      cat('\\end{verbatim}\n')
    } else {
      maxlen <- max(nchar(names(g)))
      est <- function(X, Y, x) approx(X, Y, xout=x, rule=2)$y
      z <- latexDotchart(g, names(g), auxdata=g, auxtitle='N',
                         w = 1 + est(c(2, 60), c(.5, 6), maxlen),
                         h = min(max(2.5*lg/20, 1), 8))
      cat(z, sep='\n')
    }
    cat("\n")
  }
  
  if(length(g <- object$na.detail.response)) {
    cat("\nStatistics on Response by Missing/Non-Missing Status of Predictors\n\n")
    print(unclass(g))
    cat("\n")           
  }
  if(file != '') sink()
  invisible()
}
                         
html.naprint.delete <- function(object, ...) {
  lg <- length(g <- object$nmiss)
  R <- character(0)
  if(lg) {
    if(sum(g > 0) < 4)
      R <- c('',
             'Frequencies of Missing Values Due to Each Variable<br>',
             '', '<pre>', capture.output(print(g)), '</pre>')
    else {
      maxlen <- max(nchar(names(g)))
      g  <- g[order(g)]
      fi <- tempfile(fileext='.png')
      png(fi, width=400, height=30 + length(g) * 24)
      opar <- par(mar=c(4,4,2,3), mgp=c(3-.75,1-.5,0))
      on.exit(par(opar))
      dotchart3(g, names(g), auxdata=g,
                xlab='Missing',
                main='Frequencies of NAs Due to Each Variable')
      dev.off()
      R <- c(tobase64image(fi), '<br>')
      #print(dotchartp(g, names(g), auxdata=g, auxtitle='N',
      #          main='Frequencies of Missing Values Due to Each Variable',
      #          showlegend = FALSE,
      #          sort   = 'descending',
      #          xlab   = 'Missing',
      #          width  = min(550, 300 + 20 * maxlen),
      #          height = plotlyParm$heightDotchart(lg)) ) 
    }
  }
  
  if(length(g <- object$na.detail.response)) {
    R <- c(R, '',
           'Statistics on Response by Missing/Non-Missing Status of Predictors<br>',
           '<pre>', capture.output(print(unclass(g))), '</pre>')
  }
  R
}
    
## Function to print model fit statistics
## Example:
#prStats(list('Observations', c('Log','Likelihood'),
##            c('Rank','Measures'),
##            c('Mean |difference|','Measures')),
##       list(list(N0=52, N1=48), list('max |deriv|'=1e-9,'-2 LL'=1332.23,
##            c(NA,2)),
#            list(tau=-.75, Dxy=-.64, C=.743, 2),
#            list(g=1.25, gr=11.3, 2)))
## Note that when there is an unnamed element of w, it is assumed to be
## the number of digits to the right of the decimal place (recycling of
## elements is done if fewer elements are in this vector), causing
## round(, # digits) and format(..., nsmall=# digits).  Use NA to use
## format without nsmall and without rounding (useful for integers and for
## scientific notation)

prStats <- function(labels, w, lang=c('plain', 'latex', 'html')) {
  lang  <- match.arg(lang)
  lorh  <- lang != 'plain'
  specs <- markupSpecs[[lang]]

  partial <- htmlSpecial('part')
  vbar    <- htmlTranslate('|')
  cca     <- htmlSpecial('combiningcircumflexaccent')
  beta    <- htmlGreek('beta')
  geq     <- htmlTranslate('>=')


  spaces <- function(n) if(n <= 0.5) '' else
   substring('                                                         ',
             1, floor(n))
  ## strsplit returns character(0) for ""
  ssplit <- function(x) {
    x <- strsplit(x, split='\n')
    for(i in 1 : length(x)) if(! length(x[[i]])) x[[i]] <- ''
    x
    }
  trans <- switch(lang,
                  latex = latexTranslate,
                  html  = htmlTranslate,
                  plain = function(x) x )
  ## Find maximum width used for each column
  p <- length(labels)
  width <- numeric(p)
for(i in 1:p) {
    labs <- ssplit(labels[i])[[1]]
    width[i] <- max(nchar(labs))
    u <- w[[i]]
    dig <- NA
    if(any(names(u)=='')) {
      dig <- unlist(u[names(u) == ''])
      u   <- u[names(u) != '']
    }
    lu  <- length(u)
    dig <- rep(dig, length=lu)
    fu  <- character(lu)
    for(j in 1 : length(u)) {
      uj <- u[[j]]
      nuj <- names(u)[j]
      dg <- dig[j]
      fu[j] <- if(nuj == 'Cluster on') specs$code(trans(uj))
               else
                 if(nuj == 'max |deriv|')
                   formatNP(signif(uj, 1), lang=lang)
               else
                 if(is.na(dg)) format(uj)
               else
                 if(dg < 0) formatNP(uj, -dg, pvalue=TRUE, lang=lang)
               else
                 formatNP(uj, dg, lang=lang)
    }
    names(fu) <- names(u)
    w[[i]]    <- fu
    for(j in 1 : length(u))
      width[i] <- max(width[i],
                      1 + nchar(nuj) + nchar(fu[j]))
  }
  if(lorh) {
    maxl <- max(sapply(w, length))
    z <- matrix('', nrow=maxl, ncol=p)
    fil <- if(lang == 'latex') '~\\hfill ' else htmlSpecial('emsp')

    chisq <- specs$chisq()
    
    trans <- rbind(
      'Dxy'        = c(latex = '$D_{xy}$',
                       html  = '<i>D</i><sub>xy</sub>'),
      'LR chi2'    = c(latex = paste0('LR ', chisq),
                       html  = paste0('LR ', chisq)),
      'Score chi2' = c(latex = paste0('Score ', chisq),
                       html  = paste0('Score ', chisq)),
      'Pr(> chi2)' = c(latex = 'Pr$(>\\chi^{2})$',
                       html  = paste0('Pr(', htmlTranslate('>'), chisq, ')')),
      'tau-a'      = c(latex = '$\\tau_{a}$',
                       html  = paste0(htmlGreek('tau'), '<sub>a</sub>')),
      'sigma gamma'= c(latex = '$\\sigma_{\\gamma}$',
                       html  = '&sigma;<sub>&gamma;</sub>'),
      'sigma w'    = c(latex = '$\\sigma_{w}$',
                       html  = '&sigma;<sub>w</sub>'),
      'gamma'      = c(latex = '$\\gamma$',
                       html  = htmlGreek('gamma')),
      'R2'         = c(latex = '$R^{2}$',
                       html  = '<i>R</i><sup>2</sup>'),
      'R2 adj'     = c(latex = '$R^{2}_{\\textrm{adj}}$',
                       html  = paste0('<i>R</i>', specs$subsup('adj', '2'))),
      'C'          = c(latex = '$C$',
                       html  = '<i>C</i>'),
      'g'          = c(latex = '$g$',
                       html  = '<i>g</i>'),
      'gp'         = c(latex = '$g_{p}$',
                       html  = '<i>g</i><sub>p</sub>'),
      'gr'         = c(latex = '$g_{r}$',
                       html  = '<i>g</i><sub>r</sub>'),
      'max |deriv|'   = c(latex = '$\\max|\\frac{\\partial\\log L}{\\partial \\beta}|$',
                          html  = paste0('max ', vbar, partial,
                                         'log <i>L</i>/', partial,
                                         beta, vbar)),
      'mean |Y-Yhat|' = c(latex = 'mean $|Y-\\hat{Y}|$',
                          html  = paste0('mean ', vbar, '<i>Y - Y</i>',
                                         cca, vbar)),
      'Distinct Y'   = c(latex = 'Distinct $Y$',
                       html  = 'Distinct <i>Y</i>'),
      'Median Y'   = c(latex = '$Y_{0.5}$',
                       html  = '<i>Y</i><sub>0.5</sub>'),
      '|Pr(Y>=median)-0.5|'  =
        c(latex = '$|\\overline{\\mathrm{Pr}(Y\\geq Y_{0.5})-\\frac{1}{2}}|$',
          html  = paste0('<span style="text-decoration: overline">', vbar,
                         'Pr(<i>Y</i> ', geq, ' median)-',
                       htmlSpecial('half'), vbar,
                         '</span>'))

    )
    
    for(i in 1 : p) {
      k <- names(w[[i]])
      for(j in 1 : length(k)) {
        u <- k[j]
        k[j] <- if(u %in% rownames(trans)) trans[u, lang]
        else if(grepl('R2\\(', u))   # handle R2(p,n) from R2Measures
          switch(lang,
                 plain = u,
                 latex = sub('R2\\((.*)\\)', '$R^{2}_{\\1}$', u),
                 html  = sub('R2\\((.*)\\)',
                             paste0('<i>R</i>',
                                    specs$subsup('\\1', '2')),u))
        else
          switch(lang,
                 plain = u,
                 latex = latexTranslate(u, greek=TRUE),
                 html  = htmlTranslate (u, greek=TRUE) )
      }
      z[1 : length(k), i] <- paste0(k, fil, w[[i]])
    }
    
    al <- paste0('|', paste(rep('c|', p), collapse=''))
    if(lang == 'latex')
      w <- latexTabular(z, headings=labels, align=al, halign=al,
                        translate=FALSE, hline=2, center=TRUE)
    else {
      labels <- gsub('\n', '<br>', labels)
      w <- htmlTable::htmlTable(z,
                                header=labels,
                                css.cell = 'min-width: 9em;',
                                align=al, align.header=al,
                                escape.html=FALSE)
      w <- htmltools::HTML(paste0(w, '\n'))
    }
    return(w)
  }
  z <- labs <- character(0)
  for(i in 1:p) {
    wid <- width[i]
    lab <- ssplit(labels[i])[[1]]
    for(j in 1:length(lab))
      lab[j] <- paste0(spaces((wid - nchar(lab[j])) / 2), lab[j])
    labs <- c(labs, paste(lab, collapse='\n'))
    u   <- w[[i]]
    a <- ''
    for(i in 1:length(u))
      a <- paste0(a, names(u)[i],
                 spaces(wid - nchar(u[i]) - nchar(names(u[i]))),
                 u[i],
                 if(i < length(u)) '\n')
    z <- c(z, a)
  }
  res <- rbind(labs, z)
  rownames(res) <- NULL
  capture.output(print.char.matrix(res, vsep='', hsep='    ', csep='',
                                   top.border=FALSE, left.border=FALSE))
}

## reListclean is used in conjunction with pstats
## Example:
## x <- c(a=1, b=2)
## c(A=x[1], B=x[2])
## reListclean(A=x[1], B=x[2])
## reListclean(A=x['a'], B=x['b'], C=x['c'])
## reListclean(A=x[1], B=c(x1=x[1], x2=x[2]))
## The last form causes B to be expanded into to two list elements
## named x1 and x2 and the name B is ignored
## reListclean(A=x[1], namesFrom=z) where z is only a 1 element vector will
## still override namesFrom (literally) with names(z) if 
## Update 2023-04-23: new argument dec which is appended to resulting
## vector and has elements removed if elements are removed from main
## information due to NA or NULL

#reListclean <- function(..., na.rm=TRUE) {
#  d <- list(...)
#  d <- d[sapply(d, function(x) ! is.null(x))]
#  x <- unlist(d)
#  names(x) <- names(d)
#  if(na.rm) x[! is.na(x)] else x
#}
reListclean <- function(..., dec=NULL, na.rm=TRUE) {
  d <- list(...)
  if(length(dec)) dec <- rep(dec, length=length(d))
  g <- if(na.rm) function(x) length(x) > 0 && ! all(is.na(x))
       else
         function(x) length(x) > 0
  keep <- which(sapply(d, g))
  w    <- d[keep]
  if(length(dec)) dec <- dec[keep]
  
  r <- list()
  nam <- names(w)
  i   <- 0
  nm  <- character(0)
  for(u in w) {
    i <- i + 1
    for(j in 1 : length(u)) {
      if(is.na(u[j])) next
      r <- c(r, u[j])
      nm <- c(nm, if(nam[i] != 'namesFrom' & length(u) == 1) nam[i] else {
            if(! length(names(u))) stop('vector element does not have names')
            names(u)[j] })
    }
  }
  names(r) <- nm
  c(r, dec)
}



formatNP <- function(x, digits=NULL, pvalue=FALSE,
                     lang=c('plain', 'latex', 'html')) {
  lang <- match.arg(lang)
  if(! is.numeric(x)) return(x)
  digits <- as.numeric(digits)  # Needed but can't figure out why
    x <- as.numeric(x)
  f <- if(length(digits) && ! is.na(digits))
         format(round(x, digits), nsmall=digits, scientific=1) else
         format(x, scientific=1)
  sci <- grep('e', f)
  if(length(sci)) {
    if(lang == 'latex') f[sci] <- paste0('$', latexSN(f[sci]), '$')
    else
      if(lang == 'html') f[sci] <- htmlSN(f[sci])
  }
  f <- ifelse(is.na(x), '', f)

  if(! pvalue) return(f)
  
  if(! length(digits)) stop('must specify digits if pvalue=TRUE')
  s <- ! is.na(x) & x < 10 ^ (-digits)
  if(any(s)) {
    w <- paste0('0.', paste0(rep('0', digits - 1), collapse=''), '1')
    f[s] <- switch(lang,
                   latex = paste0('\\textless ', w),
                   html  = paste0(htmlTranslate('<'), w),
                   plain = paste0('<', w))
  }
  f
}

logLik.ols <- function(object, ...) {
  ll <- getS3method('logLik', 'lm')(object)
  attr(ll, 'df') <- object$stats['d.f.'] + 2
  ll
}

logLik.rms <- function(object, ...)
  {
    dof <- unname(object$stats['d.f.'] + num.intercepts(object))
    if(inherits(object, 'psm')) dof <- dof + 1  # for sigma
    nobs <- nobs(object)
    w <- object$loglik
    if(length(w)) return(structure(w[length(w)], nobs=nobs, df=dof,
                                   class='logLik'))
    w <- object$deviance
    structure(-0.5*w[length(w)], nobs=nobs, df=dof, class='logLik')
  }

logLik.Gls <- function(object, ...) getS3method('logLik', 'gls')(object, ...)
    
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

setPb <- function(n, type=c('Monte Carlo Simulation','Bootstrap',
                       'Cross-Validation'),
                  label, usetk=TRUE, onlytk=FALSE, every=1) {
  type <- match.arg(type)
  if(!missing(label)) type <- label
  pbo <- .Options$showprogress
  if(!length(pbo)) pbo <- 'console'
  else if(is.logical(pbo)) {
    pbo <- if(pbo) 'tk' else 'none'
  }
  if(missing(every)) {
    evo <- .Options$showevery
    if(length(evo)) every <- evo
  }
  if(pbo == 'none') return(function(i, ...){invisible()})
  if(pbo == 'tk' && usetk && requireNamespace('tcltk', quietly=TRUE)) {
    pb <- tcltk::tkProgressBar(type, 'Iteration: ', 0, n)
    upb1 <- function(i, n, every, pb) {
      if(i %% every == 0)
        tcltk::setTkProgressBar(pb, i, label=sprintf('Iteration: %d', i))
      if(i == n) close(pb)
    }
    formals(upb1) <- list(i=0, n=n, every=every, pb=pb)
    return(upb1)
  }
  if(onlytk) return(function(...) {invisible()})
  upb2 <- function(i, n, every) {
    if(i %% every == 0)
      cat('Iteration: ', i, ' of ', n, '\r', sep='')
    if(i == n) cat('\n')
  }
  formals(upb2) <- list(i=0, n=n, every=every)
  upb2
}

## Function to remove one or more terms from a model formula, using
## strictly character manipulation.  This handles problems such as
## [.terms removing offset() if you subset on anything
## For each character string in which, terms like string(...) are removed.

removeFormulaTerms <- function(form, which=NULL, delete.response=FALSE) {
  if('offset' %in% which) {
    form <- formula(terms(form)[TRUE])
    which <- setdiff(which, 'offset')
  }
  ## [.terms ignores offset variables.  Above logic handles nested () unlike
  ## what is below
  form <- paste(deparse(form), collapse='')  # no string splitting
  if(delete.response) form <- gsub('.*~', '~', form)
  for(w in which) {
    pattern <- sprintf('\\+?[ ]*?%s\\(.*?\\)[ ]*?\\+{0,1}', w)  ## assume additive form
    form <- gsub(pattern, '', form)
  }
  as.formula(form)
}
