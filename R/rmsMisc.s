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

vcov.lrm <- function(object, regcoef.only=TRUE, intercepts='all', ...)
  vcov.orm(object, intercepts=intercepts, ...)

vcov.ols <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)

vcov.cph <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)

vcov.psm <- function(object, regcoef.only=TRUE, ...)
  vcov.rms(object, regcoef.only=regcoef.only, ...)

vcov.orm <- function(object, regcoef.only=TRUE, intercepts='mid', ...) {
  np   <- length(object$coefficients)
  ns   <- num.intercepts(object)
  v    <- object[['var']]
  info <- object$info.matrix
  override <- object$override_vcov_intercept
  if(length(override)) intercepts <- override
  li1  <- length(intercepts) == 1
  iat  <- if(length(v)) attr(v, 'intercepts')  # handle fit.mult.impute (?), robcov
  iref <- object$interceptRef
  i    <- c(iref, if(np > ns) (ns + 1) : np)
  if(is.numeric(intercepts) && li1 && intercepts == iref) intercepts <- 'mid'

  # Handle old fit objects from < rms 7.0-0 (only mid intercept stored) or new
  # fits run through robcov, bootcov, fit.mult.impute (all intercepts stored in var)
  if(length(v)) {
    # Old orm stored var for mid intercept only, but robcov, fit.mult.impute, bootcov
    # stored whole matrix
    type <- if(ncol(v) == np) 'full'
      else if(ncol(v) == (np - ns + 1)) 'mid'
      else stop('variance-covariance matrix stored in object$var has # intercepts ',
                'that is not 1 or all ', ns)
    if(intercepts == 'mid' && type == 'mid') return(v)
    if(intercepts == 'mid' && type == 'full') return(v[i, i, drop=FALSE])
    if(intercepts == 'all' && type == 'full') return(v)
    if(intercepts == 'none') return(if(type == 'mid') v[-1, -1, drop=FALSE]
                                    else v[-(1 : ns), -(1 : ns), drop=FALSE])
    if(is.numeric(intercepts) && type == 'full') {
      # See https://github.com/harrelfe/rms/issues/167
      i <- c(intercepts, if(np > ns) (ns + 1) : np)
      return(v[i, i, drop=FALSE])
    }
    stop('intercepts=', intercepts, ' requested for an orm fit wth $var.\n',
         'orm only stored intercept components of the variance-covariance matrix ',
         'for the middle intercept.')
  }

  # Instead of dealing with $var deal with $info

  name <- names(coef(object))
  if(is.numeric(intercepts)) i <- c(intercepts, if(np > ns) (ns + 1) : np)
  else {
    if(intercepts == 'none') i <- 'x'
    else if(intercepts == 'all') {
      v <- Matrix::as.matrix(infoMxop(info, invert=TRUE))
      dimnames(v) <- list(name, name)
      return(v)
    }
  }
  # Left with original i for mid intercept, or i just defined
  v <- Matrix::as.matrix(infoMxop(info, i=i))
  if(is.character(i) && (i == 'x')) i <- - (1 : ns)
  dimnames(v) <- list(name[i], name[i])
  v
  }

vcov.rms <- function(object, regcoef.only=TRUE, intercepts='all', ...)
  {
    np <- length(object$coefficients)
    ns <- num.intercepts(object)
    cov <- object$var
    if(length(cov)) {
      if(regcoef.only) cov <- cov[1 : np, 1 : np, drop=FALSE]
      if(length(intercepts) && intercepts == 'none' && ns > 0)
        cov <- cov[- (1 : ns), - (1 : ns), drop=FALSE]
      return(cov)
    }
    # regcoef.only does not apply to fits from lrm, orm which are the only fits using
    # info.matrix
    info <- object$info.matrix
    if(length(intercepts) && intercepts == 'none')
      infoMxop(info, i = (ns + 1) : np)
      else infoMxop(info, invert=TRUE)
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

  lr <- function(f) {
    s <- f$stats
    if(length(s)) s['Model L.R.'] else -2 * logLik(f)
    # Here Model L.R. is really shifted by null logLik but that will
    # cancel out later anyway
    # logLik is correct for glm with gaussian() family even though deviance isn't
    # glm back-calculates logLik from AIC
  }

  np <- function(f) {
    s <- f$stats
    if(length(s)) s['d.f.'] else f$rank - (any(names(coef(f)) == '(Intercept)'))
  }

  dof <- abs(np(fit2) - np(fit1))
  if(dof == 0) stop('models are not nested')

  lp1 <- length(fit1$parms);  lp2 <- length(fit2$parms)
  if(lp1 != lp2) warning('fits do not have same number of scale parameters') else
  if(lp1 == 1 && abs(fit1$parms - fit2$parms) > 1e-6)
    warning('fits do not have same values of scale parameters.\nConsider fixing the scale parameter for the reduced model to that from the larger model.')

  chisq <- abs(lr(fit2) - lr(fit1))
  r     <- c(chisq, dof, 1 - pchisq(chisq, dof))
  names(r) <- c('L.R. Chisq', 'd.f.', 'P')
  structure(list(stats    = r,
                 formula1 = formula(fit1),
                 formula2 = formula(fit2)),
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
  d <- Matrix::diag(v)^.5
  v <- Matrix::diag(Matrix::solve(v/(d %o% d)))
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

## General function to print model fit objects using latex, html, typst,
## or regular print (the default)

coef_tiny_table <- function(U, lang) {
  d <- data.frame(Term = rownames(U), unclass(U), check.names = FALSE)
  colnames(d) <- c('', colnames(U))
  nc <- ncol(d)

  x <- tinytable::tt(d, output = lang)
  x <- tinytable::format_tt(x, j = 1:nc, escape = FALSE)
  x <- tinytable::format_tt(x, i = 0,    escape = FALSE)
  ## escape=FALSE was missing entirely here -- row labels can contain
  ## Typst math (the interaction symbol, translated to $times$), and
  ## without this tinytable escapes the $ characters as literal text
  ## instead of treating them as math delimiters. stats_tiny_table
  ## already had this; coef_tiny_table simply never did.
  x <- tinytable::style_tt(x, j = 1,     align = 'l')
  x <- tinytable::style_tt(x, j = 2:nc,  align = 'r')
  x <- tinytable::style_tt(x, i = 0,     align = 'c')
  x <- tinytable::theme_striped(x)

  if(lang == 'latex')
    x <- tinytable::theme_latex(x, environment_table = FALSE)
  ## multipage=TRUE/rowhead=1 (repeated headers on long tables) were
  ## tried and dropped -- combined with environment_table=FALSE they
  ## crashed inside tinytable's internal lines_drop() (1:idx, argument
  ## of length 0). Long coefficient tables won't get automatic
  ## repeated-header page breaks the way Hmisc::latex(longtable=TRUE)
  ## did; tabularray's default tblr environment still handles overlong
  ## content, just without that feature.

  result <- tinytable::save_tt(x, lang)

  if(lang == 'latex')
    result <- paste0('\\begin{center}\n', result, '\n\\end{center}')
  else if(lang == 'typst') {
    ## Wrapping our own #block(width:100%)[#align(center)[...]] OUTSIDE
    ## the entire save_tt() result. #align(center) alone centers only
    ## within its immediate container's available width -- if that
    ## container itself has no definite width (Typst's default
    ## width:auto sizes a block to its content), there's nothing to
    ## center within. Confirmed from real output: cph/Gls/psm wrap their
    ## whole prModFit output (title, call, verbatim tables, AND both
    ## tinytable outputs) in an outer #block[ not generated by any of
    ## our own code -- likely a Quarto/knitr chunk-output artifact, not
    ## fully traced to its source -- while lrm/ols/orm do not, which is
    ## why centering worked for those and not these. Forcing our own
    ## width:100% block immediately around the content sidesteps needing
    ## to know what, if anything, surrounds it.
    typstFence <- function(x) paste0("````{=typst}\n", x, "\n````\n")
    result <- typstFence(paste0('#block(width: 100%)[#align(center)[\n',
                                result, '\n]]'))
  }
  result
}

prModFit <- function(x, title, w, digits=4, coefs=TRUE, footer=NULL,
                     lines.page=40, long=TRUE, needspace, subtitle=NULL, ...) {
  debug  <- getOption('rmsdebug', FALSE)
  if(debug) saveRDS(list(x=x, title=title, w=w, digits=digits, ...), '/tmp/prmod.rds')
  lang   <- prType()
  specs  <- markupSpecs[[lang]]
  transl <- switch(lang,
                   latex = latexTranslate,
                   html  = htmlTranslate,
                   typst = typstTranslate,
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

  ## Shared helper for verbatim/raw content under typst -- bverb()/
  ## everb() don't map cleanly here (they return separate before/after
  ## markers, but #raw() needs its content as one string argument), so
  ## call sites needing typst verbatim output use this directly instead.
  ## block: false on the #raw() element itself avoids the gray-fill
  ## show-rule Quarto's template applies to block:true raw content -- the
  ## same fix already proven for typst_verbatim() in desc.r and
  ## typst.naprint.delete's rawBlock(). But block:false is also
  ## literally inline, with no line break forced after it (same issue
  ## just fixed for catl()'s bold/center case, #strong[]/#align()[]
  ## being inline too) -- confirmed directly: content after this was
  ## running together on one line. The outer #block[...] wrapper below
  ## restores block-level separation without changing block:false on
  ## the #raw() element itself, so the gray-fill fix still holds.
  typstRaw <- function(txt) {
    txt <- gsub('\\', '\\\\', txt, fixed = TRUE)
    txt <- gsub('"', '\\"', txt, fixed = TRUE)
    raw <- paste0('#block[#raw("', txt, '", block: false)]')
    paste0("````{=typst}\n", raw, "\n````\n")
  }

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
    } else if(lang == 'typst') {
      if(verb)
        c(skipt(pre), typstRaw(paste(x, collapse = '\n')), skipt(skip))
      else if(bold || center) {
        inner <- paste(x, collapse = '\n')
        if(bold)   inner <- paste0('#strong[', inner, ']')
        if(center) inner <- paste0('#align(center)[', inner, ']')
        ## Wrapped in #block[...] -- #strong[]/#align()[] alone are
        ## inline constructs with no inherent line break after them, so
        ## whatever content followed was continuing on the same visual
        ## line regardless of blank-line separation elsewhere in the R
        ## vector. #block[] forces this to be its own distinct block.
        body <- paste0('#block[', inner, ']')
        c(skipt(pre), paste0("````{=typst}\n", body, "\n````\n"), skipt(skip))
      } else
        c(paste0(skipt(pre), x), skipt(skip))
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
    if(lang == 'typst')
      R <- c(R, typstRaw(paste(deparse(x$call), collapse = '\n')), '')
    else
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
      if(lang == 'latex') type <- 'latex.naprint.delete'
      if(lang == 'html') type <- 'html.naprint.delete'
      if(lang == 'typst') type <- 'typst.naprint.delete'
      ## tex is NOT set here for any language -- it previously was set
      ## TRUE for latex, wrapping this whole naprint.delete output
      ## (title text AND chart together) in \begin{center}...\end{center},
      ## which centered the title along with the chart -- confirmed as
      ## an unwanted side effect (title should stay left-justified,
      ## matching html/typst). Centering now happens inside
      ## latex.naprint.delete() itself, around just the chart, matching
      ## typst.naprint.delete's existing pattern.
      ##
      ## Worth watching for one side effect: with tex no longer TRUE
      ## here, `if(! tex && length(titl))` below will now emit a
      ## catl(titl, ...) title for latex if z$title happens to be set
      ## for this entry -- previously suppressed by tex=TRUE. Not
      ## something I can fully verify without seeing how the w list gets
      ## built upstream; worth confirming no unexpected extra title
      ## appears.
      ##
      ## Typst-specific note: \begin{center}/\end{center} is LaTeX-specific
      ## syntax that would break if inserted literally into a Typst
      ## document -- another reason tex must not be set TRUE for typst.
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
                 typst = x,
                 plain  = x)
        betan <- switch(lang,
                        plain = 'Beta',
                        html  = htmlGreek('beta'),
                        latex = '$\\hat{\\beta}$',
                        typst = '$hat(beta)$')

        B   <- obj$bayes
        if(length(B)) {
          U <- matrix('', nrow=nrow(B), ncol=ncol(B))
          for(i in 1:ncol(B)) {
            dig <- if(colnames(B)[i] == 'Symmetry') 2 else digits
            U[, i] <- pad(formatNP(B[, i], dig, lang=lang))
          }
          ## P(...) with P in math font, matching the rest of the table
          pn <- switch(lang, plain='P(Beta>0)',
                       html = paste0('<i>P</i>(', betan, transl('>'), '0)'),
                       latex = '$P(\\beta>0)$',
                       typst = '$P(beta>0)$')
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
                      'P(>|Z|)' =
                        formatNP(P, 4, lang=lang, pvalue=TRUE))
        if(lang == 'latex')
          colnames(U) <- c('$\\hat{\\beta}$', 'S.E.', 'Wald $Z$',
                           '$P(>|Z|)$')
        else if(lang == 'html')
            colnames(U) <- c(htmlGreek('beta'),   # did have cca
                             'S.E.', 'Wald <i>Z</i>',
                             paste0('<i>P</i>(', gt, vbar, '<i>Z</i>', vbar, ')'))
        else if(lang == 'typst')
            colnames(U) <- c('$hat(beta)$', 'S.E.', 'Wald $Z$',
                             '$P(>|Z|)$')
        if(length(errordf))
          colnames(U)[3:4] <-
            switch(lang,
                   latex = c('$t$', '$P(>|t|)$'),
                   html  = c('<i>t</i>', paste0('<i>P</i>(', gt, vbar, '<i>t</i>',
                                                vbar, ')')),
                   typst = c('$t$', '$P(>|t|)$'),
                   plain = c('t',   'P(>|t|)') )

        rownames(U) <- betanames

        if(length(obj$aux)) {
          U <- cbind(U, formatNP(obj$aux, digits, lang=lang))
          colnames(U)[ncol(U)] <- obj$auxname
        }
        }
        if(lang %in% c('latex', 'html', 'typst')) {
          R <- c(R, skipt(1))
          rownames(U) <- transl(betanames)

          if(is.numeric(coefs)) {
            U <- U[1:coefs,,drop=FALSE]
            U <- rbind(U, rep('', ncol(U)))
            rownames(U)[nrow(U)] <-
              switch(lang, html = '&hellip;', latex = '\\dots', typst = '\u2026')
          }
          ## Translate interaction symbol (*) to times symbol -- must
          ## run AFTER transl(), not before, so the rest of the term
          ## name gets properly escaped by transl() first. But the
          ## search pattern has to be format-aware: under typst,
          ## transl() (typstTranslate) has already escaped bare * to \*
          ## by this point (needed in general, since * is Typst's
          ## bold/emphasis shorthand), so the pattern must match the
          ## escaped form or the substitution leaves an orphaned
          ## backslash. Substituting before transl() instead (tried
          ## first) doesn't work either: typstTranslate() escapes EVERY
          ## $ it sees unconditionally, with no way to tell "deliberate
          ## math delimiter" from "literal text needing escaping", so a
          ## pre-inserted $times$ gets both its $ signs escaped to
          ## \$times\$. latex/html's transl() functions never touch *
          ## at all, so they still need the bare pattern.
          pat <- if(lang == 'typst') '\\*' else '*'
          rownames(U) <- gsub(pat, specs$times, rownames(U), fixed=TRUE)

          if(! missing(needspace) && lang == 'latex')
            R <- c(R, paste0('\\Needspace{', needspace, '}'))

          R <- c(R, coef_tiny_table(U, lang))
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

      w <- switch(type,
                  html.naprint.delete = do.call(type, obj),
                  typst.naprint.delete = do.call(type, obj),
                  latex.naprint.delete =
                    capture.output(do.call(type,
                                           c(obj, list(file='')))),
                   print =
                     if(lang == 'typst')
                       ## bverb()/everb() have no typst case (return
                       ## NULL), so this content was previously reaching
                       ## Pandoc completely bare -- multi-line print()
                       ## output collapsed into one run-on line
                       ## (markdown's soft-newline handling) with a
                       ## stray gray background on some lines
                       ## (Pandoc's own indented-line-as-code-block
                       ## auto-detection firing inconsistently).
                       ## typstRaw() (#raw(block:false), fenced)
                       ## sidesteps both.
                       typstRaw(paste(capture.output(do.call(type, obj)),
                                     collapse = '\n'))
                     else
                       c(bverb(),
                         capture.output(do.call(type, obj)),
                         everb()),
                   do.call(type, obj)
                   )

      R <- c(R, w,
             ## unlike do.call, eval(call(...)) dispatches on class of ...
             if(tex) '\\end{center}' else '' )
    }
  }
  if(length(footer))
    R <- c(R, paste(specs$smallskip, transl(footer)))

  if(debug)
    cat(R, sep='\n', append=TRUE, file='/tmp/rmsdebug.txt')

  switch(lang,
         html  = rendHTML(R),
         latex = rendHTML(R, html=FALSE),
         typst = rendHTML(R, html=FALSE),
         plain = cat(R, sep='\n')
         )
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
                         w = 1.25 * (1 + est(c(2, 60), c(.5, 6), maxlen)),
                         h = 1.25 * min(max(2.5*lg/20, 1), 8))
      ## Centering moved here, around just the chart -- not the title
      ## text above (cat() call a few lines up), which should stay
      ## left-justified. Previously the whole naprint.delete output was
      ## centered from outside, via prModFit's tex=TRUE dispatch flag,
      ## which centered the title along with the chart -- confirmed as
      ## an unwanted side effect.
      cat('\\begin{center}\n')
      cat(z, sep='\n')
      cat('\\end{center}\n')
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

## Typst analog of latex.naprint.delete/html.naprint.delete. Uses
## typstDotchart (native vector markup, like latexDotchart) rather than
## html.naprint.delete's PNG-embedding approach, since Typst has no need
## for a raster fallback here. Returns a plain character vector, same
## shape as html.naprint.delete -- called directly via do.call(), not
## writing to a file the way latex.naprint.delete does via sink().
typst.naprint.delete <- function(object, ...) {
  lg <- length(g <- object$nmiss)
  R <- character(0)
  typstFence <- function(x) paste0("````{=typst}\n", x, "\n````\n")
  rawBlock <- function(txt) {
    txt <- gsub('\\', '\\\\', txt, fixed = TRUE)
    txt <- gsub('"', '\\"', txt, fixed = TRUE)
    ## #block[...] wrapper: block:false on #raw() itself avoids the
    ## gray-fill show-rule, but is also literally inline with no line
    ## break forced after it -- same fix as prModFit's typstRaw().
    typstFence(paste0('#block[#raw("', txt, '", block: false)]'))
  }

  if(lg) {
    R <- c(R, 'Frequencies of Missing Values Due to Each Variable', '')
    if(sum(g > 0) < 4)
      R <- c(R, rawBlock(paste(capture.output(print(g)), collapse = '\n')))
    else {
      maxlen <- max(nchar(names(g)))
      est <- function(X, Y, x) approx(X, Y, xout = x, rule = 2)$y
      z <- typstDotchart(g, names(g), auxdata = g, auxtitle = 'N',
                         w = 1.25 * (1 + est(c(2, 60), c(.5, 6), maxlen)),
                         h = 1.25 * min(max(2.5 * lg / 20, 1), 8))
      ## Centered via the same #block(width:100%)[#align(center)[...]]
      ## pattern used for the tinytable-based tables (coef_tiny_table
      ## etc.) -- typstDotchart's own #box(...) output has no centering
      ## of its own, and align(center) alone needs a definite-width
      ## container to center within.
      R <- c(R, typstFence(paste0('#block(width: 100%)[#align(center)[\n',
                                  z, '\n]]')))
    }
  }

  if(length(g <- object$na.detail.response))
    R <- c(R, '', 'Statistics on Response by Missing/Non-Missing Status of Predictors', '',
           rawBlock(paste(capture.output(print(unclass(g))), collapse = '\n')))
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

stats_tiny_table <- function(z, labels, lang) {
  noHyphen <- function(s) {
    parts <- strsplit(s, '\n', fixed = TRUE)[[1]]
    wrapped <- vapply(parts, function(part) {
      words <- strsplit(part, ' ', fixed = TRUE)[[1]]
      w <- if(lang == 'latex') paste0('\\mbox{', words, '}')
           else if(lang == 'typst')
             paste0('#text(hyphenate: false)[', words, ']')
           else words
      paste(w, collapse = ' ')
    }, character(1))
    paste(wrapped, collapse = '\n')
  }
  labels <- vapply(labels, noHyphen, character(1))

  d <- as.data.frame(z, stringsAsFactors = FALSE)
  colnames(d) <- labels
  nc <- ncol(d)
  nr <- nrow(d)

  ## width scales with column count -- fixed values were tuned for a
  ## 3-column example and cramped a 4-column one
  ncTmp <- ncol(z)
  width <- if(lang == 'latex') min(0.95, 0.21 * ncTmp)
           else if(lang == 'typst') min(0.95, 0.233 * ncTmp)
           else NULL

  x <- if(length(width)) tinytable::tt(d, output = lang, width = width)
       else tinytable::tt(d, output = lang)

  x <- tinytable::format_tt(x, j = 1:nc, escape = FALSE)
  x <- tinytable::format_tt(x, i = 0, escape = FALSE, linebreak = "\n")
  x <- tinytable::style_tt(x, j = 1:nc, align = 'l')
  x <- tinytable::style_tt(x, i = 0, align = 'c', bold = TRUE)
  x <- tinytable::style_tt(x, i = 0, j = 1:nc, line = "b")
  x <- tinytable::style_tt(x, i = 0:nr, j = 1:nc, line = "l")
  x <- tinytable::style_tt(x, i = 0:nr, j = nc,   line = "r")
  x <- tinytable::theme_striped(x)

  if(lang == 'latex')
    x <- tinytable::theme_latex(x, environment_table = FALSE,
                                inner = "colsep=8.8pt,rowsep=1pt")
  else if(lang == 'typst')
    x <- tinytable::style_tt(x, fontsize = 0.9)

  result <- tinytable::save_tt(x, lang)

  if(lang == 'latex')
    result <- paste0('\\begin{center}\n', result, '\n\\end{center}')
  else if(lang == 'typst') {
    ## Own #block(width:100%)[#align(center)[...]] wrap outside
    ## save_tt()'s entire result -- see coef_tiny_table's identical
    ## comment for the full rationale (align(center) alone needs a
    ## definite-width container to center within; cph/Gls/psm's real
    ## output showed an outer #block[ of unknown origin, not from any of
    ## our own code, wrapping this content with no definite width).
    typstFence <- function(x) paste0("````{=typst}\n", x, "\n````\n")
    result <- typstFence(paste0('#block(width: 100%)[#align(center)[\n',
                                result, '\n]]'))
  }
  result
}

prStats <- function(labels, w, lang=c('plain', 'latex', 'html', 'typst')) {
  lang  <- match.arg(lang)
  lorh  <- lang != 'plain'
  specs <- markupSpecs[[lang]]

  partial <- htmlSpecial('part')
  vbar    <- htmlTranslate('|')
  cca     <- htmlSpecial('combiningcircumflexaccent')
  beta    <- htmlGreek('beta')
  geq     <- htmlTranslate('>=')

  debug   <- getOption('rmsdebug', FALSE)


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
                  typst = typstTranslate,
                  plain = function(x) x )
  ## Find maximum width used for each column
  if(debug) {prn(labels); prn(length(labels))}
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
    for(j in seq_len(lu)) {
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
    for(j in seq_len(length(u)))
      width[i] <- max(width[i],
                      1 + nchar(nuj) + nchar(fu[j]))
  }
  if(lorh) {
    maxl <- max(sapply(w, length))
    z <- matrix('', nrow=maxl, ncol=p)
    fil <- switch(lang,
                 latex = '~\\hfill ',
                 typst = paste0(' ', specs$hfill, ' '),
                 htmlSpecial('emsp'))

    chisq <- specs$chisq()

    trans_tbl <- rbind(
      'Dxy'        = c(latex = '$D_{xy}$',
                       html  = '<i>D</i><sub>xy</sub>',
                       typst = '$D_("xy")$'),
      'LR chi2'    = c(latex = paste0('LR ', chisq),
                       html  = paste0('LR ', chisq),
                       typst = paste0('LR ', chisq)),
      'Score chi2' = c(latex = paste0('Score ', chisq),
                       html  = paste0('Score ', chisq),
                       typst = paste0('Score ', chisq)),
      'P(> chi2)' = c(latex = '$P(>\\chi^{2})$',
                       html  = paste0('<i>P</i>(', htmlTranslate('>'), chisq, ')'),
                       typst = '$P(>chi^2)$'),
      'tau-a'      = c(latex = '$\\tau_{a}$',
                       html  = paste0(htmlGreek('tau'), '<sub>a</sub>'),
                       typst = '$tau_(a)$'),
      'sigma gamma'= c(latex = '$\\sigma_{\\gamma}$',
                       html  = '&sigma;<sub>&gamma;</sub>',
                       typst = '$sigma_(gamma)$'),
      'sigma w'    = c(latex = '$\\sigma_{w}$',
                       html  = '&sigma;<sub>w</sub>',
                       typst = '$sigma_(w)$'),
      'gamma'      = c(latex = '$\\gamma$',
                       html  = htmlGreek('gamma'),
                       typst = '$gamma$'),
      'R2'         = c(latex = '$R^{2}$',
                       html  = '<i>R</i><sup>2</sup>',
                       typst = '$R^(2)$'),
      'R2 adj'     = c(latex = '$R^{2}_{\\textrm{adj}}$',
                       html  = paste0('<i>R</i>', specs$subsup('adj', '2')),
                       typst = '$R^(2)_("adj")$'),
      'C'          = c(latex = '$C$',
                       html  = '<i>C</i>',
                       typst = '$C$'),
      'g'          = c(latex = '$g$',
                       html  = '<i>g</i>',
                       typst = '$g$'),
      'gp'         = c(latex = '$g_{p}$',
                       html  = '<i>g</i><sub>p</sub>',
                       typst = '$g_(p)$'),
      'gr'         = c(latex = '$g_{r}$',
                       html  = '<i>g</i><sub>r</sub>',
                       typst = '$g_(r)$'),
      'max |deriv|'   = c(latex = '$\\max|\\frac{\\partial\\log L}{\\partial \\beta}|$',
                          html  = paste0('max ', vbar, partial,
                                         'log <i>L</i>/', partial,
                                         beta, vbar),
                          typst = 'max $|(partial log L)/(partial beta)|$'),
      'mean |Y-Yhat|' = c(latex = 'mean $|Y-\\hat{Y}|$',
                          ## cca moved inside <i>, adjacent to the Y it
                          ## modifies -- see end-of-file notes
                          html  = paste0('mean ', vbar, '<i>Y - Y', cca,
                                         '</i>', vbar),
                          typst = 'mean $|Y-hat(Y)|$'),
      'Distinct Y'   = c(latex = 'Distinct $Y$',
                       html  = 'Distinct <i>Y</i>',
                       typst = 'Distinct $Y$'),
      'Median Y'   = c(latex = '$Y_{0.5}$',
                       html  = '<i>Y</i><sub>0.5</sub>',
                       typst = '$Y_(0.5)$'),
      '|P(Y>=median)-0.5|'  =
        c(latex = '$|\\overline{P(Y\\geq Y_{0.5})-\\frac{1}{2}}|$',
          html  = paste0('<span style="text-decoration: overline">', vbar,
                         'P(<i>Y</i> ', geq, ' median)-',
                       htmlSpecial('half'), vbar,
                         '</span>'),
          typst = '$|overline("P"(Y >= Y_(0.5)) - 1/2)|$')

    )

    for(i in seq_len(p)) {
      k <- names(w[[i]])
      for(j in seq_len(length(k))) {
        u <- k[j]
        k[j] <- if(u %in% rownames(trans_tbl)) trans_tbl[u, lang]
        else if(grepl('R2\\(', u))   # handle R2(p,n) from R2Measures
          switch(lang,
                 plain = u,
                 latex = sub('R2\\((.*)\\)', '$R^{2}_{\\1}$', u),
                 html  = sub('R2\\((.*)\\)',
                             paste0('<i>R</i>',
                                    specs$subsup('\\1', '2')),u),
                 typst = sub('R2\\((.*)\\)',
                             '$R^(2)_("\\1")$', u))
        else
          switch(lang,
                 plain = u,
                 latex = latexTranslate(u, greek=TRUE),
                 html  = htmlTranslate (u, greek=TRUE),
                 typst = typstTranslate(u, greek=TRUE) )
      }
      if(lang == 'html')
        ## HTML has no text-level equivalent to \hfill/#h(1fr) -- a
        ## plain em-space (the previous fil value) never actually
        ## right-justified the value. Needs a real flex container
        ## wrapping label and value together (open AND close tags),
        ## which doesn't fit the simple insert-a-separator-string
        ## pattern used for latex/typst -- built directly here instead
        ## of via fil.
        z[seq_len(length(k)), i] <-
          paste0('<span style="display:flex;justify-content:space-between">',
                '<span>', k, '</span><span>', w[[i]], '</span></span>')
      else
        z[seq_len(length(k)), i] <- paste0(k, fil, w[[i]])
    }

    if(lang %in% c('latex', 'html', 'typst'))
      return(stats_tiny_table(z, labels, lang))

    return(z)
  }
  z <- labs <- character(0)
  for(i in seq_len(p)) {
    wid <- width[i]
    lab <- ssplit(labels[i])[[1]]
    for(j in seq_len(length(lab)))
      lab[j] <- paste0(spaces((wid - nchar(lab[j])) / 2), lab[j])
    labs <- c(labs, paste(lab, collapse='\n'))
    u   <- w[[i]]
    a <- ''
    for(i in seq_len(length(u)))
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
    for(j in seq_len(length(u))) {
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



## typstSN: Typst analog of latexSN/htmlSN, for scientific-notation
## coefficients formatted by formatNP() below.
typstSN <- function(x) {
  x <- format(x)
  sedit(x, c("e+00", "e-0*", "e-*", "e+0*", "e+*"),
        c("", "times 10^(-*)", "times 10^(-*)", "times 10^(*)", "times 10^(*)"))
}

formatNP <- function(x, digits=NULL, pvalue=FALSE,
                     lang=c('plain', 'latex', 'html', 'typst')) {
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
    else if(lang == 'html') f[sci] <- htmlSN(f[sci])
    else if(lang == 'typst') f[sci] <- paste0('$', typstSN(f[sci]), '$')
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
                   typst = paste0(typstTranslate('<'), w),
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
## drop.terms will not remove offset()

if(FALSE) removeFormulaTerms <- function(form, which=NULL, delete.response=FALSE) {
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

## Version of removeFormulaTerms that uses the terms() and drop.terms()s Functions

removeFormulaTerms <- function(form, which=NULL, delete.response=FALSE) {
  # drop.terms will not remove offsets.  Trick it by renaming offset() terms .off.
  if('offset' %in% which) {
    form <- format(form)
    which[which == 'offset'] <- '.off.'
    z <- gsub('offset(', '.off.(', form, fixed=TRUE)
    form <- as.formula(paste(z, collapse=' '))
  }

  te <- terms(form, specials=which)
  s  <- unlist(attr(te, 'specials'))  # LHS present -> 1 added to s
  ypresent <- attr(te, 'response')
  # drop.terms counts only RHS terms
  te <- drop.terms(te, s - ypresent, keep.response= ! delete.response)
  formula(te)   # don't allow other attributes to be there
}

## -----------------------------------------------------------------------
## Consolidated notes on functions revised in this file for Typst output
## support (coef_tiny_table, prModFit, stats_tiny_table, prStats,
## typstSN, formatNP). Two different kinds of change are mixed together
## below, worth telling apart. Most of what follows is simply building
## out Typst coverage for the first time -- part of this project's goal
## from the start, not corrections to anything broken. A few items are
## genuine, unrelated bug fixes, marked as such explicitly.
##
## coef_tiny_table (new) / prModFit
## ---------------------------------
## - New helper coef_tiny_table() renders the coefficient table via
##   tinytable for lang %in% c('latex','html','typst'), replacing the
##   separate Hmisc::latex()/htmlTable::htmlTable() calls that used to
##   live inline in prModFit. Row names are an explicit first
##   (left-aligned) column; statistic columns are right-aligned; header
##   row is centered; rows are striped (theme_striped()).
## - For LaTeX, environment_table=FALSE is set explicitly on
##   theme_latex() -- not relied on as an automatic side effect of
##   multipage=TRUE, which theme_latex's own documentation claims but
##   which real output demonstrated does not actually happen (the table
##   appeared as a numbered, captioned float regardless). A manual
##   \begin{center}...\end{center} wrap is added back, since removing
##   the float also removes the \centering that came bundled with it.
##   multipage=TRUE/rowhead=1 (repeated headers on long tables) were
##   tried and dropped -- combined with environment_table=FALSE they
##   crashed inside tinytable's own internal lines_drop() (1:idx,
##   argument of length 0); long coefficient tables will not get
##   automatic repeated-header page breaks the way
##   Hmisc::latex(longtable=TRUE) did.
## - transl, pad(), and betan (implementing typst support): previously
##   handled only latex/html/plain; each gains a 'typst' case. These
##   were bare switch() calls with no default, so an unhandled language
##   would have silently returned NULL rather than erred -- worth
##   knowing for future additions of this kind.
## - Final emission switch(lang, ...) (implementing typst support):
##   gains typst = rendHTML(R, html=FALSE). Separately, a genuine bug
##   fix unrelated to typst: the latex case was previously
##   cat(R, sep='\n') rather than rendHTML(R, html=FALSE) -- this was
##   the actual root cause of needing results='asis' in chunks that
##   mixed prModFit-based printing with other latex() calls; corrected
##   here alongside the typst addition.
## - Pr(...) notation changed to P(...), with P set in math font,
##   consistently across latex/html/typst (was plain-text "Pr" followed
##   by math-mode parenthetical content in the latex/html cases).
## - Not yet extended to typst this pass, flagged rather than left
##   silently incomplete: bverb()/everb() and catl()'s verb=TRUE path
##   fall through to plain-style (unstyled) output under typst -- no
##   monospace box around the deparsed call/printed objects.
##
## typst.naprint.delete (new)
## ---------------------------------
## - New function, Typst analog of latex.naprint.delete/
##   html.naprint.delete, added to prModFit's naprint.delete dispatch
##   (type <- 'typst.naprint.delete' when lang=='typst', tex stays FALSE
##   -- \begin{center}/\end{center}, what tex=TRUE triggers, is
##   LaTeX-specific and would break in a Typst document) and to the
##   switch(type, ...) call site (do.call(type, obj) directly, same
##   shape as html.naprint.delete). Uses typstDotchart for the
##   >=4-variable missing-data chart -- native vector markup, like
##   latexDotchart, rather than html.naprint.delete's PNG-embedding
##   approach. Verbatim print() output (the <4-variable case, and the
##   na.detail.response section) uses #raw(block:false), fenced, the
##   same pattern established for prModFit's 'print' type branch.
##
## stats_tiny_table (new) / prStats
## ---------------------------------
## - New helper stats_tiny_table() renders prStats' final table via
##   tinytable for lang %in% c('latex','html','typst'), replacing
##   latexTabular()/htmlTable::htmlTable(). This also removed a stray
##   Helvetica-font-wrapping artifact that had been latexTabular()'s
##   default behavior.
## - match.arg(lang), the trans switch, and two generic
##   label-translation switches (implementing typst support): each
##   gains a 'typst' case, previously latex/html/plain only.
## - Every row of the internal trans_tbl translation matrix
##   (implementing typst support) gains a typst entry -- previously
##   latex/html only, and since trans_tbl[u, lang] indexes by column
##   name, an unhandled row would have errored outright rather than
##   returned something usable. FIRST DRAFT, not individually
##   compile-tested the way most of this effort has been -- one entry
##   has since been confirmed and corrected by testing: 'max |deriv|'
##   originally used 'diff' as the Typst symbol name for the
##   partial-derivative sign, which does not exist; 'partial' is
##   correct and is now in place. The remaining highest-risk entry is
##   the overline-wrapped median-deviation row (a long expression inside
##   overline(), not yet exercised).
## - fil (implementing typst support, confirmed via dedicated compile
##   testing at this exact multi-column table shape): gains a typst
##   case using specs$hfill ('#h(1fr)'), the Typst analog of LaTeX's
##   \hfill, requiring format_tt(escape=FALSE) and an explicit tt(width=)
##   to have any visible effect -- confirmed empirically, not assumed.
## - P(> chi2) row specifically: Pr(...) -> P(...), math font, across
##   latex/html/typst, matching the same change in prModFit. The
##   trans_tbl row key itself was also corrected from 'Pr(> chi2)' to
##   'P(> chi2)' -- the actual runtime stat name cph.s builds (renamed
##   from 'Score P' before calling prStats) has no "r", so the original
##   key (present unchanged since the first uploaded version of this
##   file) never actually matched, silently falling through to the
##   generic typstTranslate(u, greek=TRUE)/equivalent path for every
##   language, not just typst -- confirmed by comparing real output
##   against a working row ('Score chi2') using the same lookup
##   mechanism.
## - Several rounds of styling refinement, applied via tinytable's
##   documented style_tt()/theme_*() API rather than hand-built
##   per-format markup: header row centered, bold, with a bottom rule
##   (style_tt(i=0, ...)); data rows left-aligned (each cell already
##   handles its own internal left/right split via fil/hfill); vertical
##   column borders via style_tt(line="l"/"r") (theme_grid(), tried
##   first, added unwanted interior row lines instead of just a header
##   rule); row striping via theme_striped(); reduced row spacing via
##   tabularray's rowsep option for LaTeX, and an approximate fontsize
##   reduction for Typst (no confirmed direct row-spacing parameter
##   exists there); no-hyphenation in headers via per-WORD
##   \mbox{}/#text(hyphenate:false) wrapping (wrapping whole multi-word
##   segments was tried first, but that blocks ordinary word-wrap too,
##   which caused long headers to overflow into neighboring columns
##   instead of wrapping to a second line); table width scaled by
##   column count rather than a single fixed value, since the fixed
##   value was tuned only for a 3-column example and cramped a 4-column
##   one.
## - 'mean |Y-Yhat|' row, html entry: a genuine bug, not a typst gap --
##   the combining circumflex accent (cca) was placed after the closing
##   </i> tag rather than immediately adjacent to the second Y it's
##   meant to modify. Unicode combining characters only combine with a
##   preceding base character when directly adjacent, with no
##   intervening markup -- across an HTML tag boundary, a browser
##   renders the accent as its own standalone glyph instead, which is
##   what a lone combining circumflex accent looks like: a caret to the
##   right ("Y^") rather than a hat over the Y ("Ŷ"). Confirmed directly
##   from real rendered output (latex/typst were already correct, only
##   html showed the caret). Fixed by moving cca inside the <i> tag,
##   directly adjacent to the Y.
##
## typstSN (new) / formatNP
## ---------------------------------
## - New helper typstSN(), the Typst analog of latexSN()/htmlSN(), for
##   scientific-notation coefficients.
## - match.arg(lang) and the scientific-notation branch (implementing
##   typst support): gain a 'typst' case, using typstSN().
## - pvalue branch (implementing typst support): gains a typst case
##   using typstTranslate('<'), mirroring how the html case already
##   uses htmlTranslate('<') -- needed because Typst reads a bare '<' as
##   the start of label-reference syntax.
## -----------------------------------------------------------------------
