Gls <-
  function (model, data = sys.frame(sys.parent()), correlation = NULL, 
    weights = NULL, subset, method = c("REML", "ML"), na.action = na.omit, 
    control = list(), verbose = FALSE, B=0, dupCluster=FALSE,
            pr=FALSE, opmeth=c('optimize','optim'), x=FALSE)
            
{
    require(nlme)
    glsEstimate <- nlme:::glsEstimate
    
    Call <- match.call()
    opmeth <- match.arg(opmeth)
    controlvals <- glsControl()
    if (!missing(control))
      {
        if (!is.null(control$nlmStepMax) && control$nlmStepMax < 0) {
          warning("Negative control$nlmStepMax - using default value")
          control$nlmStepMax <- NULL
        }
        controlvals[names(control)] <- control
      }
    if (!inherits(model, "formula") || length(model) != 3)
        stop("\nModel must be a formula of the form \"resp ~ pred\"")

    method <- match.arg(method)
    REML <- method == "REML"
    if (length(correlation))
      groups <- getGroupsFormula(correlation)
    else groups <- NULL
    glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
    mfArgs <- list(formula = asOneFormula(formula(glsSt), model, 
                     groups), data = data, na.action = na.action)
    if (!missing(subset)) 
      mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]

    mfArgs$drop.unused.levels <- TRUE
    dataMod <- do.call("model.frame", mfArgs)
    rn <- origOrder <- row.names(dataMod)  ## rn FEH 6apr03
    if (length(groups))
      {
        groups <- eval(parse(text = paste("~1", deparse(groups[[2]]), 
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
    if(!length(dul) || dul)
      {
        on.exit(options(drop.unused.levels=dul))
        options(drop.unused.levels=FALSE)
      }
    X <- Design(X)
    atrx <- attributes(X)
    desatr <- atrx$Design
    mt <- atrx$terms
    attr(X,'Design') <- NULL
    
    contr <- lapply(X, function(el) if (inherits(el, "factor")) 
                    contrasts(el))
    contr <- contr[!unlist(lapply(contr, is.null))]
    X <- model.matrix(model, X)
    dimnames(X)[[2]] <- cn <- c('Intercept',desatr$colnames)
    y <- eval(model[[2]], dataMod)
    N <- nrow(X)
    p <- ncol(X)
    parAssign <- attr(X, "assign")
    fTerms <- terms(as.formula(model))
    namTerms <- attr(fTerms, "term.labels")
    if (attr(fTerms, "intercept") > 0)
      namTerms <- c("Intercept", namTerms)
    namTerms <- factor(parAssign, labels = namTerms)
    parAssign <- split(order(parAssign), namTerms)
    ## Start FEH 4apr03
    if(B > 0)
      {
        bootcoef <- matrix(NA, nrow=B, ncol=p, dimnames=list(NULL,cn))
        bootcorr <- numeric(B)
        Nboot    <- integer(B)
        if(length(grps))
          {
            obsno    <- split(1:N,grps)
            levg     <- levels(grps)
            ng       <- length(levg)
            if(!length(levg)) stop('program logic error')
          }
        else
          {
            obsno <- 1:N
            levg  <- NULL
            ng    <- N
          }
      }
    for(j in 0:B)
      {
        if(j == 0) s <- 1:N else
        {
          if(ng == N) s <- sample(1:N, N, replace=TRUE)
          else
            {
              grps.sampled <- sample(levg, ng, replace=TRUE)
              s <- unlist(obsno[grps.sampled])
              dataMods <- dataMod[s,]
              if(!dupCluster)
                {
                  grp.freqs <- table(grps)
                  newgrps <- factor(rep(paste('C',1:ng,sep=''),
                                        table(grps)[grps.sampled]))
                  dataMods$id <- newgrps
                }
            }
          Nboot[j] <- Nb <- length(s)
          if(pr) cat(j,'\r')
        }
        attr(glsSt, "conLin") <-
          if(j==0)
            list(Xy = array(c(X, y), c(N, p + 1),
                   list(rn, c(cn, deparse(model[[2]])))), 
                 dims = list(N = N, p = p, REML = as.integer(REML)),
                 logLik = 0)
          else
            list(Xy = array(c(X[s,,drop=FALSE], y[s]), c(Nb, p + 1),
                   list(rn[s], c(cn, deparse(model[[2]])))), 
                 dims = list(N = Nb, p = p, REML = as.integer(REML)),
                 logLik = 0)
        ## FEH colnames(X) -> cn, ncol(X) -> p, j>0 case above
        glsEstControl <- controlvals[c("singular.ok", "qrTol")]
        glsSt <- Initialize(glsSt, if(j==0) dataMod else dataMods,
                            glsEstControl)
        parMap <- attr(glsSt, "pmap")
        numIter <- numIter0 <- 0
        repeat
          {
            oldPars <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
            if (length(coef(glsSt)))
              {
                co <- c(coef(glsSt))  ## FEH
                if(opmeth == 'optimize' && co > 1) opmeth <- 'optim'
                best <- if(opmeth=='optimize')
                  optimize(function(z)
                           -logLik(glsSt,z), lower=-12, upper=12)$minimum else
                optim(fn = function(glsPars)
                      -logLik(glsSt, glsPars),
                      par = co,  ## FEH
                      method = "BFGS", 
                      control = list(trace = controlvals$msVerbose, 
                        reltol = if (numIter == 0)
                        controlvals$msTol else 100 * 
                        .Machine$double.eps,
                        maxit = controlvals$msMaxIter))$par
                coef(glsSt) <- best   ## FEH
              }
            attr(glsSt, "glsFit") <- glsEstimate(glsSt, control = glsEstControl)
            if (!needUpdate(glsSt)) 
              break
            numIter <- numIter + 1
            glsSt <- update(glsSt, if(j==0) dataMod else dataMods)  ## FEH
            aConv <- c(attr(glsSt, "glsFit")[["beta"]], coef(glsSt))
            conv <- abs((oldPars - aConv)/ifelse(aConv == 0, 1, aConv))
            aConv <- c(beta = max(conv[1:p]))
            conv <- conv[-(1:p)]
            for (i in names(glsSt))
              {
                if (any(parMap[, i]))
                  {
                    aConv <- c(aConv, max(conv[parMap[, i]]))
                    names(aConv)[length(aConv)] <- i
                  }
              }
            if (verbose)
              {
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
      if(j > 0)
        {
          bootcoef[j,] <- attr(glsSt, "glsFit")[["beta"]]
          bootcorr[j]  <- coef(glsSt$corStruct, unconstrained=FALSE)
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
    if (length(grps))
      {
        grps <- grps[revOrder]
        Fitted <- Fitted[revOrder]
        Resid <- y[revOrder] - Fitted
        attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt)[revOrder])
      }
    else
      {
        Resid <- y - Fitted
        attr(Resid, "std") <- glsFit$sigma/(varWeights(glsSt))
      }
    attr(Resid, 'label') <- 'Residuals'
    cr <- class(Resid)
    if(length(cr) && any(cr == 'labelled'))
      {
        if(length(cr) == 1) Resid <- unclass(Resid) else
        class(Resid) <- setdiff(cr, 'labelled')
      }
    
    if (controlvals$apVar && FALSE)
      apVar <-
        nlme:::glsApVar(glsSt, glsFit$sigma,
                        .relStep  = controlvals[[".relStep"]],  
                        minAbsPar = controlvals[["minAbsParApVar"]],
                        natural   = controlvals[["natural"]])
    else
      apVar <- "Approximate variance-covariance matrix not available"
    dims <- attr(glsSt, "conLin")[["dims"]]
    dims[["p"]] <- p
    attr(glsSt, "conLin") <- NULL
    attr(glsSt, "glsFit") <- NULL
    estOut <- list(modelStruct = glsSt, dims = dims, contrasts = contr, 
                   coefficients = glsFit[["beta"]], varBeta = varBeta,
                   sigma = glsFit$sigma, g=GiniMd(Fitted), 
                   apVar = apVar, logLik = glsFit$logLik,
                   numIter = if(needUpdate(glsSt)) numIter else numIter0,  
                   groups = grps, call = Call, method = method,
                   fitted = Fitted,  
                   residuals = Resid, parAssign = parAssign,
                   Design=desatr, assign=DesignAssign(desatr, 1, mt),
                   formula=model, terms=fTerms, opmeth=opmeth,
                   B=B, boot.Coef=if(B > 0) bootcoef,
                   boot.Corr=if(B > 0) bootcorr,
                   Nboot=if(B > 0) Nboot,
                   var=if(B > 0) var(bootcoef),
                   x=if(x) X)
    
    ## Last 2 lines FEH 29mar03
    if (inherits(data, "groupedData"))
      {
        attr(estOut, "units") <- attr(data, "units")
        attr(estOut, "labels") <- attr(data, "labels")
      }
    attr(estOut, "namBetaFull") <- colnames(X)
    estOut$fitFunction <- 'Gls'
    class(estOut) <- c('Gls','rms','gls')
    estOut
  }


print.Gls <- function(x, digits=4, coefs=TRUE, latex=FALSE, ...)
{
  ## Following taken from print.gls with changes marked FEH

  summary.gls <- nlme:::summary.gls

  k <- 0
  z <- list()
  
  dd <- x$dims
  errordf <- dd$N - dd$p
  mCall <- x$call
  Title <- if (inherits(x, "gnls")) "Generalized Nonlinear Least Squares Fit"
  else paste("Generalized Least Squares Fit by",
             ifelse(x$method == "REML", "REML", "Maximum Mikelihood"))

  ltype <- if (inherits(x, "gnls")) 'Log-likelihood' else
   paste('Log-', ifelse(x$method == "REML", "restricted-", ""),
         'likelihood',  sep='')
  if(latex) ltype <- paste(ltype, ' ', sep='')
  
  misc <- reVector(Obs=dd$N,
                   Clusters=if(length(x$groups)) length(unique(x$groups)) else
                    dd$N,
                   g=x$g)
  
  llike <- reVector(ll=x$logLik,
                    'Model d.f.' = dd$p - 1,
                    sigma  = x$sigma,
                    'd.f.' = errordf)
  names(llike)[1] <- ltype
  k <- k + 1
  z[[k]] <- list(type='stats',
                 list(
                      headings = list('', ''),
                      data     = list(c(misc, c(NA,NA,3)),
                                      c(llike, c(2,NA,digits,NA)))))

  if(any(names(x)=='var') && length(x$var))
    {
      se <- sqrt(diag(x$var))
      beta <- coef(x)
      k <- k + 1
      z[[k]] <- list(type='coefmatrix',
                     list(coef = beta, se= se),
                     title='Using bootstrap variance estimates')
    }
  else
    {
      s <- summary.gls(x)$tTable
      k <- k + 1
      z[[k]] <- list(type='coefmatrix',
                     list(coef = s[,'Value'], se = s[,'Std.Error'],
                          errordf = errordf))
    }
  

  if (length(x$modelStruct) > 0)
    {
      k <- k + 1
      z[[k]] <- list(type='print', list(summary(x$modelStruct)))
    }

  if(x$B > 0)
    {
      k <- k + 1
      z[[k]] <- list(type='cat', list('Bootstrap repetitions:',x$B))

      tn <- table(x$Nboot)
      if(length(tn) > 1)
        {
          k < k + 1
          z[[k]] <- list(type='print', list(tn), 
                         title = 'Table of Sample Sizes used in Bootstraps')
        }
      else
        {
          k <- k + 1
          z[[k]] <- list(type='cat',
                         list('Bootstraps were all balanced with respect to clusters'))
        }
      
      dr <- diag(x$varBeta)/diag(x$var)
      k <- k + 1
      z[[k]] <- list(type='print', list(round(dr, 2)),
                     title = 'Ratio of Original Variances to Bootstrap Variances')
      k <- k + 1
      r <- round(quantile(x$boot.Corr, c(.025,.975)),3)
      names(r) <- c('Lower','Upper')
      z[[k]] <- list(type='print', list(r),
                     title = 'Bootstrap Nonparametric 0.95 Confidence Limits for Correlation Parameter')
    }
  
  prModFit(x, title=Title, z, digits=digits, coefs=coefs, latex=latex, ...)
  invisible()
}

vcov.Gls <- function(object, ...)
  if(any(names(object)=='var') && length(object$var))
  object$var else object$varBeta

predict.Gls <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","cterms", "ccterms", "adjto",
             "adjto.data.frame", "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=type=="terms", ...)
  {
    type <- match.arg(type)
    predictrms(object, newdata, type, se.fit, conf.int, conf.type,
               incl.non.slopes, non.slopes, kint,
               na.action, expand.na, center.terms, ...)
  }
    

latex.Gls <- function(...) latexrms(...)
