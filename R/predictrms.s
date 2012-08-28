##newdata=data frame, vector,  matrix, or list.  All but first assume data
##need coding, e.g. categorical variables are given as integers
##variable missing for all obs -> use adjust-to value in limits
##(means (parms) for matrx)

predictrms <-
  function(fit, newdata=NULL,
           type=c("lp","x","data.frame","terms","cterms","ccterms","adjto",
             "adjto.data.frame","model.frame"),
           se.fit=FALSE, conf.int=FALSE,
           conf.type=c('mean','individual','simultaneous'),
           incl.non.slopes=NULL, non.slopes=NULL, kint=1,
           na.action=na.keep, expand.na=TRUE,
           center.terms=type=="terms", ref.zero=FALSE, ...)
{
  type <- match.arg(type)
  conf.type <- match.arg(conf.type)
  if(conf.type == 'simultaneous') {
    require(multcomp)
    if(missing(newdata) || !length(newdata))
      stop('newdata must be given if conf.type="simultaneous"')
  }
  ## R does not preserve missing(x)
  mnon.slopes <- missing(non.slopes) || !length(non.slopes)  # was missing( )

  at        <- fit$Design
  assume    <- at$assume.code
  Limval    <- Getlim(at, allow.null=TRUE, need.all=FALSE)
  Values    <- Limval$values
  non.ia    <- assume != 9
  non.strat <- assume != 8
  f <- sum(non.ia)
  nstrata   <- sum(assume==8)
  somex     <- any(non.strat)
  rnam      <- NULL
  cox <- inherits(fit, "cph") ||
          (length(fit$fitFunction) && any(fit$fitFunction=='cph'))
  naa <- fit$na.action
  if(!expand.na)
    naresid <- function(a,b) b #don't really call naresid if drop NAs

  parms <- at$parms
  name  <- at$name
  coeff <- fit$coefficients
  nrp   <- num.intercepts(fit)

  if(mnon.slopes)
    {
      non.slopes <- rep(0,nrp)
      non.slopes[kint] <- 1
    }
  else if(nrp>0 & length(non.slopes)!=nrp)
    stop("length of non.slopes incompatible with fit")

  int.pres <- nrp > 0
  if(somex) cov <- vcov(fit, regcoef.only=TRUE)    #remove scale params
  if(missing(incl.non.slopes) || !length(incl.non.slopes))
    incl.non.slopes <- !mnon.slopes | (!missing(kint)) | 
                        int.pres | type!="x"
  int.pres <- int.pres & incl.non.slopes

  assign <- fit$assign

  nama <- names(assign)[1]
  asso <- 1*(nama=="Intercept" | nama=="(Intercept)")

  Center <- if(cox)fit$center else 0

  oldopts <- options(contrasts=c(factor="contr.treatment",ordered="contr.poly"),
                     Design.attr=at)

  ## In SV4 options(two lists) causes problems
  on.exit({options(contrasts=oldopts$contrasts)
           options(Design.attr=NULL)})

  Terms <- delete.response(terms(formula(fit), specials='strat'))
  
  ## Remove any offset terms - in rms they have been omitted from newdata
  off <- attr(Terms, 'offset')
  if(length(off)) Terms <- drop.terms(Terms, off)
  
  attr(Terms,"response")  <- 0
  attr(Terms,"intercept") <- 1
  ##Need intercept whenever design matrix is generated to get
  ##current list of dummy variables for factor variables
  stra <- attr(Terms, "specials")$strat

  Terms.ns <- if(length(stra)) Terms[-stra] else Terms

  if(conf.int)
    {
      vconstant <- 0
      if(conf.type=='individual')
        {
          vconstant <- fit$stats['Sigma']^2
          if(is.na(vconstant))
            stop('conf.type="individual" requires that fit be from ols')
        }
      zcrit <- if(length(idf <- fit$df.residual)) qt((1+conf.int)/2, idf) else
      qnorm((1+conf.int)/2)
    }

  if(type != "adjto" & type != "adjto.data.frame")
    {
      X <- NULL
      if(missing(newdata) || !length(newdata))
        {
          if(type=="lp" && length(fit$linear.predictors))
            {
              LP <- naresid(naa, fit$linear.predictors)
              if(kint > 1)
                LP <- LP - fit$coefficients[1] + fit$coefficients[kint]
              if(length(stra <- fit$Strata))
                attr(LP, "strata") <- naresid(naa, stra)
              if(!se.fit && !conf.int)return(LP)
              else
                if(length(fit$se.fit))
                  {
                    if(kint>1)
                      warning("se.fit is retrieved from the fit but it corresponded to kint=1")
                    retlist <- list(linear.predictors=LP)
                    if(se.fit) retlist$se.fit <-
                      naresid(naa, fit$se.fit)
                    if(conf.int)
                      {
                        plminus <- zcrit*sqrt(retlist$se.fit^2 + vconstant)
                        retlist$lower <- LP - plminus
                        retlist$upper <- LP + plminus
                      }
                    return(retlist)
                  }
            }
          else
            if(type=="x") return(structure(naresid(naa, fit$x),
                 strata=if(length(stra <- fit$Strata))
                         naresid(naa, stra) else NULL))
          X <- fit$x
          rnam <- dimnames(X)[[1]]
          if(!any(names(fit)=="x")) X <- NULL  #fit$x can get fit$xref
          if(!length(X))
            stop("newdata not given and fit did not store x")
        }
      if(!length(X))
        {
          if(!is.data.frame(newdata))
            {
              if(is.list(newdata))
                {
                  loc <- if(!length(names(newdata))) 1:f else name[assume!=9]
                  new <- matrix(double(1),
                                nrow=length(newdata[[1]]),
                                ncol=length(newdata))
                  for(j in 1:ncol(new)) new[,j] <- newdata[[loc[j]]]
                  newdata <- new
                }
              if(!is.matrix(newdata)) newdata <- matrix(newdata, ncol=f)
              if(ncol(newdata) != f)
                stop("# columns in newdata not= # factors in design")
              X  <- list()
              k  <- 0
              ii <- 0
              for(i in (1:length(assume))[non.ia])
                {
                  ii <- ii+1
                  xi <- newdata[,ii]
                  as <- assume[i]
                  allna <- all(is.na(xi))
                  if(as==5 | as==8)
                    {
                      xi <- as.integer(xi)
                      levels(xi) <- parms[[name[i]]]
                      class(xi) <- "factor"
                    }
                  else if(as==10)
                    {
                      if(i==1) ifact <- 1
                      else ifact <- 1 + sum(assume[1:(i-1)]!=8)
                      ##	Accounts for assign not being output for strata factors
                      ncols <- length(assign[[ifact+asso]])
                      if(allna)
                        {
                          xi <- matrix(double(1),
                                       nrow=length(xi), ncol=ncols)
                          for(j in 1:ncol(xi)) xi[,j] <- parms[[name[i]]][j]
                        }
                      else xi <- matrix(xi, nrow=length(xi), ncol=ncols)
                    }
                  ##	Duplicate single value for all parts of matrix
                  k <- k + 1
                  X[[k]] <- xi
                }
              names(X) <- name[non.ia]
              attr(X, "row.names") <- as.character(1:nrow(newdata))
              class(X) <- "data.frame"
              newdata <- X
              ##Note: data.frame() converts matrix variables to individual variables
              if(type=="data.frame") return(newdata)
            }
          else
            {
              ## Need to convert any factors to have all levels in original fit
              ## Otherwise, wrong dummy variables will be generated by model.matrix
              nm <- names(newdata)
              for(i in 1:ncol(newdata))
                {
                  j <- match(nm[i], name)
                  if(!is.na(j))
                    {
                      asj <- assume[j]
                      w <- newdata[,i]
                      V <- NULL
                      if(asj==5 | asj==7 | asj==8 | 
                         (name[j] %in% names(Values) &&
                          length(V <- Values[[name[j]]]) && is.character(V)))
                        {
                          if(length(Pa <- parms[[name[j]]])) V <- Pa
                          newdata[,i] <- factor(w, V)
                          ## Handles user specifying numeric values without quotes, that
                          ## are levels
                          ww <- is.na(newdata[,i]) & !is.na(unclass(w))
                          if(any(ww)) 	{
                            cat("Error in predictrms: Values in",names(newdata)[i],
                                "not in",V,":\n")
                            print(as.character(w[ww]),quote=FALSE); stop()
                          }
                        }
                    }
                }
            }
          X <- model.frame(Terms, newdata, na.action=na.action, ...)
          if(type=="model.frame") return(X)
          naa <- attr(X, "na.action")
          rnam <- row.names(X)
          
          strata <- list()
          nst <- 0
          ii <- 0
          for(i in 1:ncol(X)) {
            ii <- ii + 1
            xi <- X[[i]]
            asi <- attr(xi, "assume.code")
            as <- assume[ii]
            if(!length(asi) && as==7) {
              attr(X[,i],"contrasts") <- 
                attr(scored(xi,name=name[ii]),"contrasts")
              if(length(xi)==1) warning("a bug in model.matrix can produce incorrect results\nwhen only one observation is being predicted for an ordered variable")
            }
            
              if(as==8) {
                nst <- nst+1
                ste <- paste(name[ii], parms[[name[ii]]], sep='=')
                strata[[nst]] <- factor(ste[X[,i]], ste)
              }
          }
          X <- if(!somex) NULL
          else
            if(int.pres && nrp==1) model.matrix(Terms.ns, X)
            else
              model.matrix(Terms.ns, X)[,-1,drop=FALSE]
          
          if(nstrata > 0) {
            names(strata) <- paste("S",1:nstrata,sep="")
            strata <- interaction(as.data.frame(strata), drop=FALSE)
          }
        }
      else strata <- attr(X,"strata")
      added.col <- if(incl.non.slopes & (nrp>1 | !int.pres)) nrp else 0
      if(incl.non.slopes & nrp > 0 & somex & added.col > 0) {
        xx <- matrix(double(1),
                     nrow=nrow(X), ncol=added.col)
        for(j in 1:nrp) xx[,j] <- non.slopes[j]
      }
      else xx <- NULL
    }
  
  ## For models with multiple intercepts, delete elements of covariance matrix
  ## containing unused intercepts
  elements.to.delete <- 9999
  if(somex && nrp>1) {
    i <- (1:nrp)[non.slopes==0]; cov <- cov[-i,-i,drop=FALSE] 
    elements.to.delete <- i
  }
  
  if(type=="adjto" | type=="adjto.data.frame" | ref.zero |
     (center.terms && type %in% c("terms","cterms","ccterms")) | 
     (cox & (se.fit | conf.int))) {
    ## Form design matrix for adjust-to values
    adjto <- list()
    ii <- 0
    for(i in (1:length(assume))[non.ia]) {
      ii <- ii+1
      xi <- Getlimi(name[i], Limval, need.all=TRUE)[2]
      if(assume[i]==5 | assume[i]==8)
        xi <- factor(xi,parms[[name[i]]])
      else
        if(assume[i]==7) xi <- scored(xi, name=name[i])
        else
          if(assume[i]==10)
            xi <- matrix(parms[[name[i]]],nrow=1) #matrx col medians
      adjto[[ii]] <- xi
    }
    names(adjto) <- name[non.ia]
    attr(adjto,"row.names") <- "1"
    class(adjto) <- "data.frame"
    if(type=="adjto.data.frame") return(adjto)
    adjto <- model.frame(Terms, adjto)
    adjto <- if(int.pres) model.matrix(Terms.ns, adjto) else
    model.matrix(Terms.ns,adjto)[, -1, drop=FALSE]
    
    if(type=="adjto") {
      k <- if(int.pres) 1:length(coeff) else (nrp+1):length(coeff)
      if(is.matrix(adjto))
        dimnames(adjto) <- list(dimnames(adjto)[[1]],names(coeff)[k])
      else names(adjto) <- names(coeff)[k]
      return(adjto)
    }
  }
  
  if(length(xx) && type %nin% c("terms","cterms","ccterms") &&
     incl.non.slopes) {
    X <- cbind(xx, X)
    dimnames(X) <- list(rnam, names(coeff))
    if(cox & (se.fit | conf.int)) adjto <- c(xx[1,], adjto)
  }
  
  else if(somex) dimnames(X) <- 
    list(rnam,names(coeff)[(1 + length(coeff) - ncol(X)):length(coeff)])

  if(type=="x") return(
       structure(naresid(naa,X),
                 strata=if(nstrata > 0)  naresid(naa,strata) else NULL,
                 na.action=if(expand.na)NULL else naa)
       )
  
  if(type=="lp") {
    if(somex) {
      if(any(elements.to.delete==9999)) cof <- coeff
      else {
        cof <- coeff[-elements.to.delete]
        X <- X[,-elements.to.delete,drop=FALSE]
      }
      xb <- matxv(X, cof) - Center
      names(xb) <- rnam
    }
    else {
      xb <- numeric(0)
      if(nstrata > 0) attr(xb, 'strata') <- naresid(naa, strata)
      return(structure(if(se.fit) list(linear.predictors=xb,
                                       se.fit=rep(NA, length(xb))) else
                       xb,
                       na.action=if(expand.na) NULL else naa))
    }
    xb <- naresid(naa, xb)
    if(nstrata > 0) attr(xb,"strata") <- naresid(naa,strata)
    ycenter <- if(ref.zero && somex) matxv(adjto, cof) - Center else 0
    
    if(ref.zero || ((se.fit || conf.int) && somex)) {
      if(cox || ref.zero) X <- sweep(X, 2, adjto) #Center columns
      se <- drop(sqrt(((X %*% cov) * X) %*% rep(1, ncol(X))))
      names(se) <- rnam
      
      sef <- naresid(naa, se)
      ww <- if(conf.int || se.fit) {
        if(se.fit)
          list(linear.predictors = xb - ycenter, se.fit = sef) else
        list(linear.predictors = xb - ycenter)
      }
      else xb - ycenter
      retlist <- structure(ww, 
                           na.action=if(expand.na) NULL else naa)
      if(conf.int) {
        if(conf.type == 'simultaneous') {
          u <- confint(glht(fit, X,
                            df=if(length(idf)) idf else 0),
                       level=conf.int)$confint
          retlist$lower <- u[,'lwr']
          retlist$upper <- u[,'upr']
        } else {
          plminus <- zcrit*sqrt(sef^2 + vconstant)
          retlist$lower <- xb - plminus - ycenter
          retlist$upper <- xb + plminus - ycenter
        }
      }
      return(retlist)
    }
    else return(structure(xb - ycenter, na.action=if(expand.na)NULL else naa))
  }

  if(type %in% c("terms", "cterms", "ccterms")) {
    if(!somex)
      stop('type="terms" may not be given unless covariables present')
    
    usevar <- if(type=="terms") non.strat else rep(TRUE, length(assume))
    fitted <- array(0, c(nrow(X), sum(usevar)),
                    list(rnam, name[usevar]))
    if(se.fit) se <- fitted
    if(center.terms) {
      if(ncol(adjto) != ncol(X)) {
        if(dimnames(adjto)[[2]][1] %in% c('Intercept','(Intercept)') &&
           dimnames(X)[[2]][1]    %nin% c('Intercept','(Intercept)'))
          adjto <- adjto[,-1,drop=FALSE]
        if(ncol(adjto) != ncol(X)) stop('program logic error')
      }
      X <- sweep(X, 2, adjto) # center columns
    }
    num.intercepts.not.in.X <- length(coeff) - ncol(X)
    j <- 0
    for(i in (1:length(assume))[usevar]) {
      j <- j+1
      if(assume[i]!=8) { # non-strat factor; otherwise leave fitted=0
        k <- assign[[j+asso]]
        ko <- k - num.intercepts.not.in.X
        fitted[,j] <- matxv(X[,ko,drop=FALSE], coeff[k])
        if(se.fit) se[,j] <-
          (((X[, ko, drop=FALSE]  %*% cov[ko, ko, drop=FALSE]) * 
            X[, ko, drop=FALSE]) %*% rep(1, length(ko)))^.5
      }
    }
      if(type=="cterms") {
        ## Combine all related interation terms with main effect terms
        w <- fitted[, non.ia, drop=FALSE]  # non-interaction terms
        for(i in 1:f) {
          ia <- interactions.containing(at, i) # subscripts of interaction terms related to predictor i
          if(length(ia)) w[, i] <- rowSums(fitted[, c(i,ia), drop=FALSE])
        }
        fitted <- w
      }

      if(type=='ccterms') {
        z <- combineRelatedPredictors(at)
        f <- length(z$names)
        w <- matrix(NA, ncol=f, nrow=nrow(fitted))
        colnames(w) <- sapply(z$names, paste, collapse=', ')
        for(i in 1:f)
          w[,i] <- rowSums(fitted[,z$namesia[[i]], drop=FALSE])
        fitted <- w
      }
    
    fitted <- structure(naresid(naa, fitted),
                        strata=if(nstrata==0) NULL else naresid(naa, strata))
    
  if(se.fit) {
    return(structure(list(fitted=fitted, se.fit=naresid(naa,se)),
                     na.action=if(expand.na)NULL else naa)) 	}
  else return(structure(fitted, na.action=if(expand.na)NULL else naa))
  }
}   
