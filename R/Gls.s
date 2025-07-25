## This is a modification of the gls function in the nlme package.
## gls is authored by Jose Pinheiro, Douglas Bates, Saikat DebRoy,
## Deepayan Sarkar, and R-core

Gls <-
  function (model, data = sys.frame(sys.parent()), correlation = NULL, 
    weights = NULL, subset, method = c("REML", "ML"), na.action = na.omit, 
    control = list(), verbose = FALSE, B=0, dupCluster=FALSE,
            pr=FALSE, x=FALSE) {
    
    Call <- match.call()
    controlvals <- glsControl()
    if (!missing(control)) {
      if (!is.null(control$nlmStepMax) && control$nlmStepMax < 0) {
        warning("Negative control$nlmStepMax - using default value")
        control$nlmStepMax <- NULL
      }
      controlvals[names(control)] <- control
    }
    if (!inherits(model, "formula") || length(model) != 3L)
      stop("\nModel must be a formula of the form \"resp ~ pred\"")

    method <- match.arg(method)
    REML <- method == "REML"
    if (! is.null(correlation))
      groups <- getGroupsFormula(correlation)
    else groups <- NULL
    glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
    model <- terms(model, data=data)  ## new
    mfArgs <- list(formula = asOneFormula(formula(glsSt), model, 
                     groups), data = data, na.action = na.action)
    if (!missing(subset)) 
      mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2L]]
    
    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call("model.frame", mfArgs)
    rn <- origOrder <- row.names(dataMod)  ## rn FEH 6apr03
    if (length(groups)) {
      groups <- eval(parse(text = paste("~1", deparse(groups[[2L]]), 
                                        sep = "|")))
      grps <- getGroups(dataMod, groups,
                        level = length(getGroupsFormula(groups, 
                                                        asList = TRUE)))
      ord <- order(grps)
      grps <- grps[ord]
      dataMod <- dataMod[ord, , drop = FALSE]
      rn <- rn[ord]
      revOrder <- match(origOrder, rn)
    }
    else grps <- NULL
    X <- model.frame(model, dataMod)
    dul <- .Options$drop.unused.levels
    if(!length(dul) || dul) {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }
    X <- Design(X)
    atrx <- attributes(X)
    sformula <- atrx$sformula
    desatr <- atrx$Design
    mt <- atrx$terms
    mmcolnames <- desatr$mmcolnames
    
    attr(X,'Design') <- NULL
    
    contr <- lapply(X, function(el) if (inherits(el, "factor")) 
                                      contrasts(el))
    contr <- contr[!unlist(lapply(contr, is.null))]
    X <- model.matrix(model, X)
    parAssign <- attr(X, "assign")
    fixedSigma <- controlvals$sigma > 0
    alt <- attr(mmcolnames, 'alt')
    if(! all(mmcolnames %in% colnames(X)) && length(alt)) mmcolnames <- alt
    X <- X[, c('(Intercept)', mmcolnames), drop=FALSE]
    colnames(X) <- cn <- c('Intercept', desatr$colnames)
    
    y <- eval(model[[2L]], dataMod)
    N <- nrow(X)
    p <- ncol(X)
    fTerms <- terms(as.formula(model))
    namTerms <- attr(fTerms, "term.labels")
    if (attr(fTerms, "intercept") > 0) namTerms <- c("Intercept", namTerms)
    namTerms <- factor(parAssign, labels = namTerms)
    parAssign <- split(order(parAssign), namTerms)
    ## Start FEH 4apr03
    if(B > 0) {
      bootcoef <- matrix(NA, nrow=B, ncol=p, dimnames=list(NULL,cn))
      Nboot    <- integer(B)
      if(length(grps)) {
        obsno    <- split(1L : N,grps)
        levg     <- levels(grps)
        ng       <- length(levg)
        if(!length(levg)) stop('program logic error')
      } else {
        obsno <- 1L : N
        levg  <- NULL
        ng    <- N
      }
    }
    for(j in 0 : B) {
      if(j == 0) s <- 1L : N else {
                               if(ng == N) {
                                 s <- sample(1L : N, N, replace=TRUE)
                                 dataMods <- dataMod[s,]
                               }
          else {
            grps.sampled <- sample(levg, ng, replace=TRUE)
            s <- unlist(obsno[grps.sampled])
            dataMods <- dataMod[s,]
            if(!dupCluster) {
              grp.freqs <- table(grps)
              newgrps <- factor(rep(paste('C',1L : ng,sep=''),
                                    table(grps)[grps.sampled]))
              dataMods$id <- newgrps
            }
          }
                               Nboot[j] <- Nb <- length(s)
                               if(pr) cat(j,'\r')
                             }
      attr(glsSt, "conLin") <-
        if(j==0)
          list(Xy = array(c(X, y), c(N, p + 1L),
                          list(rn, c(cn, deparse(model[[2L]])))), 
               dims = list(N = N, p = p, REML = as.integer(REML)),
               logLik = 0, sigma=controlvals$sigma,
               fixedSigma=fixedSigma)
        else
          list(Xy = array(c(X[s,,drop=FALSE], y[s]), c(Nb, p + 1L),
                          list(rn[s], c(cn, deparse(model[[2L]])))), 
               dims = list(N = Nb, p = p, REML = as.integer(REML)),
               logLik = 0, sigma=controlvals$sigma,
               fixedSigma=fixedSigma)
      ## FEH colnames(X) -> cn, ncol(X) -> p, j>0 case above
      glsEstControl <- controlvals[c("singular.ok", "qrTol")]
      ## qrTol above not in gls
      glsSt <- Initialize(glsSt, if(j==0) dataMod else dataMods,
                          glsEstControl)
      parMap <- attr(glsSt, "pmap")
      numIter <- numIter0 <- 0
      repeat {
        co <- c(coef(glsSt))  ## FEH
        oldPars <- c(attr(glsSt, "glsFit")[["beta"]], co)
        if (length(co)) {
          optRes <- if(controlvals$opt == 'nlminb') {
            nlminb(co, function(glsPars) -logLik(glsSt, glsPars),
                   control = list(trace = controlvals$msVerbose, 
                                  iter.max = controlvals$msMaxIter))
          }
          else {
            optim(co, function(glsPars) -logLik(glsSt, glsPars),
                  method = controlvals$optimMethod, 
                  control = list(trace = controlvals$msVerbose, 
                                 maxit = controlvals$msMaxIter,
                                 reltol = if (numIter ==  0)
                                            controlvals$msTol else 100 * .Machine$double.eps))
          }
          coef(glsSt) <- optRes$par
        }
        else optRes <- list(convergence = 0)
        attr(glsSt, "glsFit") <- glsEstimate(glsSt, control = glsEstControl)
        if (!needUpdate(glsSt)) {
          if (optRes$convergence) stop(optRes$message)
          break
        }
        numIter <- numIter + 1L
        glsSt <- update(glsSt, if(j==0) dataMod else dataMods)  ## FEH
        aConv <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
        conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1L, aConv))
        aConv <- c(beta = max(conv[1L : p]))
        conv <- conv[-(1L : p)]
        for (i in names(glsSt)) {
          if (any(parMap[, i])) {
            aConv <- c(aConv, max(conv[parMap[, i]]))
            names(aConv)[length(aConv)] <- i
          }
        }
        if (verbose) {
          cat("\nIteration:", numIter)
          ## cat("\nObjective:", format(aNlm$value), "\n")
          ## ERROR: aNlm doesn't exist.  Need to fix.
          print(glsSt)
          cat("\nConvergence:\n")
          print(aConv)
        }
        if (max(aConv) <= controlvals$tolerance) break
        if (numIter > controlvals$maxIter)
          stop("Maximum number of iterations reached without convergence.")
      }
      if(j > 0) {
        bootcoef[j,]  <- attr(glsSt, "glsFit")[["beta"]]
        bootc <- coef(glsSt$corStruct, unconstrained=FALSE)
        if(j == 1L) {
          ncb <- ncol(bootc)
          if(!length(ncb)) ncb <- length(bootc)
          bootcorr <- matrix(NA, nrow=B, ncol=ncb,
                             dimnames=list(NULL, names(bootc)))
        }
        bootcorr[j,] <- bootc
      }
      if(j==0) glsSt0 <- glsSt  ## FEH 4apr03
    }  ## end bootstrap reps
    if(pr && B > 0) cat('\n')
    glsSt <- glsSt0   ## FEH
    
    glsFit <- attr(glsSt, "glsFit")
    namBeta <- names(glsFit$beta)
    attr(parAssign, "varBetaFact") <- varBeta <- glsFit$sigma * 
      glsFit$varBeta * sqrt((N - REML * p)/(N - p))
    varBeta <- crossprod(varBeta)
    dimnames(varBeta) <- list(namBeta, namBeta)
    Fitted <- fitted(glsSt)
    if (length(grps)) {
      grps <- grps[revOrder]
      Fitted <- Fitted[revOrder]
      Resid <- y[revOrder] - Fitted
      attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt)[revOrder])
    } else {
      Resid <- y - Fitted
      attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt))
    }
    names(Resid) <- names(Fitted) <- origOrder
    attr(Resid, 'label') <- 'Residuals'
    cr <- class(Resid)
    if(length(cr) && any(cr == 'labelled')) {
      if(length(cr) == 1L) Resid <- unclass(Resid) else
                                                     class(Resid) <- setdiff(cr, 'labelled')
    }
    
    # fixedSimga should be set before `glsApVar`
    attr(glsSt, 'fixedSigma') <- fixedSigma
    apVar <- if (isTRUE(controlvals$apVar))    # see https://github.com/harrelfe/rms/issues/157
               glsApVar(glsSt, glsFit$sigma,
                        .relStep  = controlvals[[".relStep"]],
                        minAbsPar = controlvals[["minAbsParApVar"]],
                        natural   = controlvals[["natural"]])
             else
               "Approximate variance-covariance matrix not available"
    dims <- attr(glsSt, "conLin")[["dims"]]
    dims[["p"]] <- p
    attr(glsSt, "conLin") <- NULL
    attr(glsSt, "glsFit") <- NULL
    grpDta <- inherits(data, 'groupedData')
    estOut <- list(modelStruct = glsSt, dims = dims, contrasts = contr, 
                   coefficients = glsFit[["beta"]], varBeta = varBeta,
                   sigma = if(fixedSigma) controlvals$sigma else glsFit$sigma,
                   g=GiniMd(Fitted), 
                   apVar = apVar, logLik = glsFit$logLik,
                   numIter = if(needUpdate(glsSt)) numIter else numIter0,  
                   groups = grps, call = Call, method = method,
                   fitted = Fitted,  
                   residuals = Resid, parAssign = parAssign,
                   Design=desatr, assign=DesignAssign(desatr, 1L, mt),
                   formula=model, sformula=sformula, terms=fTerms,
                   B=B, boot.Coef=if(B > 0) bootcoef,
                   boot.Corr=if(B > 0) bootcorr,
                   Nboot=if(B > 0) Nboot,
                   var=if(B > 0) var(bootcoef),
                   x=if(x) X[, -1L, drop=FALSE])
    
    ## Last 2 lines FEH 29mar03
    if(grpDta) {
      attr(estOut, "units")  <- attr(data, "units")
      attr(estOut, "labels") <- attr(data, "labels")
    }
    attr(estOut, "namBetaFull") <- colnames(X)
    class(estOut) <- c('Gls','rms','gls')
    estOut
  }


print.Gls <- function(x, digits=4, coefs=TRUE,
                      title, ...)
{
  ## Following taken from print.gls with changes marked FEH

  summary.gls <- getS3method('summary', 'gls')

  k <- 0
  z <- list()
  
  dd <- x$dims
  errordf <- dd$N - dd$p
  mCall <- x$call
  if(missing(title))
    title <- if (inherits(x, "gnls")) "Generalized Nonlinear Least Squares Fit"
    else paste("Generalized Least Squares Fit by",
               ifelse(x$method == "REML", "REML", "Maximum Mikelihood"))

  ltype <- if (inherits(x, "gnls")) 'Log-likelihood' else
   paste('Log-', ifelse(x$method == "REML", "restricted-", ""),
         'likelihood',  sep='')
  if(prType() == 'latex') ltype <- paste(ltype, ' ', sep='')
  
  misc <- reListclean(Obs=dd$N,
                   Clusters=if(length(x$groups)) length(unique(x$groups)) else
                    dd$N,
                   g=x$g,
                   dec=c(NA, NA, 3L))
  
  llike <- reListclean(ll=x$logLik,
                    'Model d.f.' = dd$p - 1L,
                    sigma  = x$sigma,
                    'd.f.' = errordf,
                    dec=c(2L,NA,digits,NA))
  names(llike)[1L] <- ltype
  k <- k + 1L
  z[[k]] <- list(type='stats',
                 list(
                      headings = c('', ''),
                      data     = list(misc, llike)))

  if(any(names(x)=='var') && length(x$var))
    {
      se <- sqrt(diag(x$var))
      beta <- coef(x)
      k <- k + 1L
      z[[k]] <- list(type='coefmatrix',
                     list(coef = beta, se= se),
                     title='Using bootstrap variance estimates')
    }
  else
    {
      ## summary.gls calls BIC which tries to use logLik.rms.
      ## Make it use logLik.gls instead
      class(x) <- 'gls'
      s <- summary.gls(x)$tTable
      k <- k + 1L
      z[[k]] <- list(type='coefmatrix',
                     list(coef = s[,'Value'], se = s[,'Std.Error'],
                          errordf = errordf))
    }
  

  if (length(x$modelStruct) > 0)
    {
      k <- k + 1L
      z[[k]] <- list(type='print', list(summary(x$modelStruct)))
    }

  if(x$B > 0)
    {
      k <- k + 1L
      z[[k]] <- list(type='cat', list('Bootstrap repetitions:',x$B))

      tn <- table(x$Nboot)
      if(length(tn) > 1L)
        {
          k < k + 1L
          z[[k]] <- list(type='print', list(tn), 
                         title = 'Table of Sample Sizes used in Bootstraps')
        }
      else
        {
          k <- k + 1L
          z[[k]] <- list(type='cat',
                         list('Bootstraps were all balanced with respect to clusters'))
        }
      
      dr <- diag(x$varBeta)/diag(x$var)
      k <- k + 1L
      z[[k]] <- list(type='print', list(round(dr, 2L)),
                     title = 'Ratio of Original Variances to Bootstrap Variances')
      k <- k + 1L
      r <- round(t(apply(x$boot.Corr, 2L, quantile, probs=c(.025,.975))), 3L)
      colnames(r) <- c('Lower','Upper')
      z[[k]] <- list(type='print', list(r),
                     title = 'Bootstrap Nonparametric 0.95 Confidence Limits for Correlation Parameters')
    }
  
  prModFit(x, title=title, z, digits=digits, coefs=coefs, ...)
}

vcov.Gls <- function(object, intercepts='all', ...)
  {
    v <- if(any(names(object)=='var') && length(object$var))
      object$var else object$varBeta
    if(length(intercepts) == 1L && intercepts == 'none')
      v <- v[-1L, -1L, drop=FALSE]
    v
  }

predict.Gls <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","cterms", "ccterms", "adjto",
             "adjto.data.frame", "model.frame"),
           se.fit=FALSE, conf.int=FALSE,
           conf.type=c('mean','individual','simultaneous'),
           kint=1L,
           na.action=na.keep, expand.na=TRUE, center.terms=type=="terms", ...)
  {
    type <- match.arg(type)
    predictrms(object, newdata, type, se.fit, conf.int, conf.type,
               kint=kint,
               na.action, expand.na, center.terms, ...)
  }
    

latex.Gls <- function(..., file='', inline=FALSE, append=FALSE) {
  z <- latexrms(..., inline=inline)
  if(inline) return(z)
  if(file == '' && prType() != 'plain') return(rendHTML(z, html=FALSE))
  cat(z, file=file, append=append, sep='\n')
  invisible()
}
