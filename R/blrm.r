##' Bayesian Binary and Ordinal Logistic Regression
##'
##' Uses \code{rstan} with pre-compiled Stan code whose location is given by the user in \code{options(stancompiled='...')} to get posterior draws of parameters from a binary logistic or proportional odds semiparametric ordinal logistic model.  The Stan code internally using the qr decompositon on the design matrix so that highly collinear columns of the matrix do not hinder the posterior sampling.  The parameters are transformed back to the original scale before returning results to R.   Design matrix columns re centered before running Stan, so Stan diagnostic output will have the intercept terms shifted but the results of \code{blrm} for intercepts are for the original uncentered data.  The only prior distributions for regression betas are normal with mean zero, and the vector of prior standard deviations is given in \code{priorsd}.  These priors are for the qr-projected design matrix elements, except that the very last element is not changed.  So if one has a single non-interactive linear or binary variable for which a skeptical prior is designed, put that variable last in the model.
##'
##' The partial proportional odds model of Peterson and Harrell (1990) is implemented, and is invoked when the user specifies a second model formula as the \code{ppo} argument.  This formula has no left-hand-side variable, and has right-side variables that are a subset of those in \code{formula} specifying for which predictors the proportional odds assumption is relaxed.
##' 
##' \code{blrm} also handles single-level hierarchical random effects models for the case when there are repeated measurements per subject which are reflected as random intercepts, and a different experimental model that allows for AR(1) serial correlation within subject.  For both setups, a \code{cluster} term in the model signals the existence of subject-specific random effects, and an additional model term \code{aTime(time variable)} signals the use of the AR(1) within-subject model.  The \code{time} variable must be integer valued and there can be arbitrary gaps between measurements.  However if the maximum time exceeds 200 or so one can expect much longer computation time.  When \code{aTime()} is present, the cluster-specific random effects then become the random effect for a subject's \code{time=1} record.  When \code{aTime} is specified, the covariates must not change over records within subject (no time-dependent covariates), as only the covariate vector for the earliest time for a subject is used.
##'
##' See \url{https://hbiostat.org/R/rms/blrm.html} for multiple examples with results.
##' @title blrm
##' @param formula a R formula object that can use \code{rms} package enhancements such as the restricted interaction operator
##' @param ppo formula specifying the model predictors for which proportional odds is not assumed
##' @param keepsep a single character string containing a regular expression applied to design matrix column names, specifying which columns are not to be QR-orthonormalized, so that priors for those columns apply to the original parameters.  This is useful for treatment and treatment interaction terms.  For example \code{keepsep='treat'} will keep separate all design matrix columns containing \code{'treat'} in their names.  Some characters such as the caret used in polynomial regression terms will need to be escaped by a double backslash.
##' @param data a data frame
##' @param subset a logical vector or integer subscript vector specifying which subset of data whould be used
##' @param na.action default is \code{na.delete} to remove missings and report on them
##' @param priorsd vector of prior standard deviations.  If the vector is shorter than the number of model parameters, it will be repeated until the length equals the number of parametertimes.
##' @param priorsdppo vector of prior standard deviations for non-proportional odds parameters.  As with \code{priorsd} the last element is the only one for which the SD corresponds to the original data scale.
##' @param conc the Dirichlet distribution concentration parameter for the prior distribution of cell probabilities at covariate means.  The default is the reciprocal of 0.8 + 0.35 max(k, 3) where k is the number of Y categories.  The default is chosen to make the posterior mean of the intercepts more closely match the MLE.  For optimizing, the concentration parameter is always 1.0 to obtain results very close to the MLE for providing the posterior mode.
##' @param psigma defaults to 1 for a half-t distribution with 4 d.f., location parameter \code{rsdmean} and scale parameter \code{rsdsd}
##' @param rsdmean the assumed mean of the prior distribution of the standard deviation of random effects.  When \code{psigma=2} this is the mean of an exponential distribution and defaults to 1.  When \code{psigma=1} this is the mean of the half-t distribution and defaults to zero.
##' @param rsdsd applies only to \code{psigma=1} and is the scale parameter for the half t distribution for the SD of random effects, defaulting to 1.
##' @param ar1sdmean the assumed mean of the prior distribution of the standard deviation of within-subject white noise.   The setup is the same as with \code{rsdmean}.
##' @param iter number of posterior samples per chain for [rstan::sampling] to run
##' @param chains number of separate chains to run
##' @param refresh see [rstan::sampling].  The default is 0, indicating that no progress notes are output.  If \code{refresh > 0} and \code{progress} is not \code{''}, progress output will be appended to file \code{progress}.  The default file name is \code{'stan-progress.txt'}.
##' @param progress see \code{refresh}.  Defaults to \code{''} if \code{refresh = 0}.  Note: If running interactively but not under RStudio, \code{rstan} will open a browser window for monitoring progress.
##' @param x set to \code{FALSE} to not store the design matrix in the fit.  \code{x=TRUE} is needed if running \code{blrmStats} for example.
##' @param y set to \code{FALSE} to not store the response variable in the fit
##' @param loo set to \code{FALSE} to not run \code{loo} and store its result as object \code{loo} in the returned object.  \code{loo} defaults to \code{FALSE} if the sample size is greater than 1000, as \code{loo} requires the per-observation likelihood components, which creates a matrix N times the number of posterior draws.
##' @param ppairs set to a file name to run \code{rstan::pairs} and store the resulting png plot there.  Set to \code{TRUE} instead to directly plot these diagnostics.  The default is not to run \code{pairs}.
##' @param method set to \code{'optimizing'} to run the Stan optimizer and not do posterior sampling, \code{'both'} (the default) to run both the optimizer and posterior sampling, or \code{'sampling'} to run only the posterior sampling and not compute posterior modes. Running \code{optimizing} is a way to obtain maximum likelihood estimates and allows one to quickly study the effect of changing the prior distributions.  When \code{method='optimizing'} is used the result returned is not a standard \code{blrm} object but is instead the parameter estimates, -2 log likelihood, and optionally the Hession matrix (if you specify \code{hessian=TRUE} in ...).  When \code{method='both'} is used, \code{rstan::sampling} and \code{rstan::optimizing} are both run, and parameter estimates (posterior modes) from \code{optimizing} are stored in a matrix \code{param} in the fit object, which also contains the posterior means and medians, and other results from \code{optimizing} are stored in object \code{opt} in the \code{blrm} fit object.  When random effects are present, \code{method} is automatically set to \code{'sampling'} as maximum likelihood estimates without marginalizing over the random effects do not make sense.
##' @param inito intial value for optimization.  The default is the \code{rstan} default \code{'random'}.  Frequently specifying \code{init=0} will benefit when the number of distinct Y categories grows or when using \code{ppo} hence 0 is the default for that.
##' @param inits initial value for sampling, defaults to \code{inito}
##' @param standata set to \code{TRUE} to return the Stan data list and not run the model
##' @param debug set to \code{TRUE} to output timing and progress information to /tmp/debug.txt
##' @param file a file name for a \code{saveRDS}-created file containing or to contain the saved fit object.  If \code{file} is specified and the file does not exist, it will be created right before the fit object is returned, less the large \code{rstan} object.  If the file already exists, its stored \code{md5} hash string \code{datahash} fit object component is retrieved and compared to that of the current \code{rstan} inputs.  If the data to be sent to \code{rstan}, the priors, and all sampling and optimization options and stan code are identical, the previously stored fit object is immediately returned and no new calculatons are done.
##' @param ... passed to \code{rstan:optimizing}.  The \code{seed} parameter is a popular example.
##' @return an \code{rms} fit object of class \code{blrm}, \code{rmsb}, \code{rms} that also contains \code{rstan} results under the name \code{rstan}.  In the \code{rstan} results, which are also used to produce diagnostics, the intercepts are shifted because of the centering of columns of the design matrix done by \code{blrm}.  With \code{method='optimizing'} a class-less list is return with these elements: \code{coefficients} (MLEs), \code{beta} (non-intercept parameters on the QR decomposition scale), \code{deviance} (-2 log likelihood), \code{return_code} (see \code{rstan::optimizing}), and, if you specified \code{hessian=TRUE} to \code{blrm}, the Hessian matrix.  To learn about the scaling of orthogonalized QR design matrix columns, look at the \code{xqrsd} object in the returned object.  This is the vector of SDs for all the columns of the transformed matrix.  Those kept out by the \code{keepsep} argument will have their original SDs.
##' @examples
##' \dontrun{
##'   options(stancompiled='~/R/stan')    # need this always
##'   stanCompile()    # do this once per computer to compile centrally
##'   getHdata(Titanic3)
##'   dd <- datadist(titanic3); options(datadist='dd')
##'   f <- blrm(survived ~ (rcs(age, 5) + sex + pclass)^2, data=titanic3)
##'   f                   # model summary using print.blrm
##'   coef(f)             # compute posterior mean parameter values
##'   coef(f, 'median')   # compute posterior median values
##'   stanDx(f)           # print basic Stan diagnostics
##'   s <- stanGet(f)     # extract rstan object from fit
##'   plot(s, pars=f$betas)       # Stan posteriors for beta parameters
##'   stanDxplot(s)       # Stan diagnostic plots by chain
##'   blrmStats(f)        # more details about predictive accuracy measures
##'   ggplot(Predict(...))   # standard rms output
##'   summary(f, ...)     # invokes summary.rms
##'   contrast(f, ...)    # contrast.rms computes HPD intervals
##'   plot(nomogram(f, ...)) # plot nomogram using posterior mean parameters
##'
##'   # Fit a random effects model to handle multiple observations per
##'   # subject ID
##'   f <- blrm(outcome ~ rcs(age, 5) + sex + cluster(id), data=mydata)
##'
##'   # Fit a random effects model that respects serial correlation within
##'   # subject ID
##'   f <- blrm(outcome ~ rcs(age, 5) + sex + cluster(id) + time(visit))
##' } 
##' @author Frank Harrell and Ben Goodrich
##' @seealso \code{\link{print.blrm}}, \code{\link{blrmStats}}, \code{\link{stanDx}}, \code{\link{stanGet}}, \code{\link{coef.rmsb}}, \code{\link{vcov.rmsb}}, \code{\link{print.rmsb}}, \code{\link{coef.rmsb}}, [stanCompile]
##' @md
blrm <- function(formula, ppo=NULL, keepsep=NULL,
                 data, subset, na.action=na.delete,
								 priorsd=rep(100, p), priorsdppo=rep(100, pppo),
                 conc=1./(0.8 + 0.35 * max(k, 3)),
                 psigma=1, rsdmean=if(psigma == 1) 0 else 1,
                 rsdsd=1, ar1sdmean=1,
								 iter=2000, chains=4, refresh=0,
                 progress=if(refresh > 0) 'stan-progress.txt' else '',
								 x=TRUE, y=TRUE, loo=n <= 1000, ppairs=NULL,
                 method=c('both', 'sampling', 'optimizing'),
                 inito=if(length(ppo)) 0 else 'random', inits=inito,
                 standata=FALSE, file=NULL, debug=FALSE,
                 ...) {

  if(missing(data))   data <- environment(formula)
  msubset <- missing(subset)
  
	call <- match.call()
  method <- match.arg(method)

  prevhash <- NULL
  if(length(file) && file.exists(file)) {
    prevfit  <- readRDS(file)
    prevhash <- prevfit$datahash
    }

  if(debug) debug <- function(...)
    cat(..., format(Sys.time()), '\n', file='/tmp/debug.txt', append=TRUE)
  else debug <- function(...) {}

  ## Function to drop only the ith dimension of a 3-array when subscripting
  ## that dimension to the jth element where j is a single integer
  ## result is a matrix
  ## Note: a[,,j, drop=TRUE] drops dimensions other than the 3rd
  ## when they have only one element
  dropadim <- function(a, i, j) {
    n <- dimnames(a)
    d <- dim(a)
    b <- if(i == 1)      a[j,  ,  , drop=FALSE]
         else if(i == 2) a[ , j,  , drop=FALSE]
         else            a[,   , j, drop=FALSE]
    d <- d[-i]
    dn <- if(length(n)) n[-i]
    matrix(as.vector(b), nrow=d[1], ncol=d[2], dimnames=dn)
    }
  
  nact <- NULL

  tform   <- terms(formula, specials=c('cluster', 'aTime'), data=data)
  yname   <- as.character(formula[2])
  
  dul <- .Options$drop.unused.levels
  if(!length(dul) || dul) {
    on.exit(options(drop.unused.levels=dul))
    options(drop.unused.levels=FALSE)
  }

  en <- environment(formula)
  assign(envir=en, 'aTime', function(x) x)
  
  requireNamespace('rstan', quietly=TRUE)

  m <- if(msubset) model.frame(formula,
                               data       = data,
                               na.action = na.action,
                               drop.unused.levels=TRUE)
       else
         model.frame(formula,
                     data=data,
                     subset=subset,
                     na.action=na.action,
                     drop.unused.levels=TRUE)

  omit <- attr(m, 'na.action')$omit

  X <- Design(m)   # Design handles cluster()
  
  cluster     <- attr(X, 'cluster')
  clustername <- attr(X, 'clustername')
  time        <- attr(X, 'time')
  timename    <- attr(X, 'timename')

  if(length(time) & ! length(cluster))
    stop('may not specify aTime() without cluster()')
  
  atrx       <- attributes(X)
  sformula   <- atrx$sformula
  nact       <- atrx$na.action
  Terms      <- atrx$terms
  attr(Terms, "formula") <- formula
  atr        <- atrx$Design
  mmcolnames <- atr$mmcolnames

  Y <- model.extract(X, 'response')
  offs <- atrx$offset
  if(!length(offs)) offs <- 0
  X <- model.matrix(Terms, X)
  alt <- attr(mmcolnames, 'alt')
  if(! all(mmcolnames %in% colnames(X)) && length(alt)) mmcolnames <- alt
  X <- X[, mmcolnames, drop=FALSE]
  colnames(X) <- atr$colnames
  notransX <- if(length(keepsep)) grep(keepsep, atr$colnames)
  if(length(keepsep) & ! length(notransX))
    warning('keepsep did not apply to any column of X')
  notransXn <- atr$colnames[notransX]

  Z <- notransZn <- NULL

  if(length(ppo)) {
    if(length(omit)) data <- data[- omit, ]
    m <- if(msubset) model.frame(ppo,
                                 data      = data,
                                 na.action = na.fail,
                                 drop.unused.levels=TRUE)
         else
           model.frame(ppo,
                       data=data,
                       subset=subset,
                       na.action=na.fail,
                       drop.unused.levels=TRUE)

    Z <- Design(m)
  
    zatrx       <- attributes(Z)
    zsformula   <- atrx$sformula
    zTerms      <- zatrx$terms
    attr(zTerms, "formula") <- ppo
    zatr        <- zatrx$Design
    mmcolnames  <- zatr$mmcolnames

    Z <- model.matrix(zTerms, Z)
    alt <- attr(mmcolnames, 'alt')
    if(! all(mmcolnames %in% colnames(Z)) && length(alt)) mmcolnames <- alt
    Z <- Z[, mmcolnames, drop=FALSE]
    colnames(Z) <- zatr$colnames
    notransZ  <- if(length(keepsep)) grep(keepsep, zatr$colnames)
    notransZn <- zatr$colnames[notransZ]
}

  if(! length(time)) {
    wqrX  <- selectedQr(X, center=TRUE, not=notransX)
    Xs    <- wqrX$X
    xbar  <- wqrX$xbar
    xqrsd <- apply(Xs, 2, sd)
    wqrX$X <- NULL

    if(length(ppo)) {
      wqrZ   <- selectedQr(Z, center=TRUE, not=notransZ)
      Zs     <- wqrZ$X
      zbar   <- wqrZ$xbar
      wqrZ$X <- NULL
      }
    }
	
	n    <- nrow(X)
	p    <- ncol(X)

  ynumeric <- is.numeric(Y)
  if(ynumeric) {
    mediany <- quantile(Y, probs=.5, type=1L)
    yu <- sort(unique(Y))
    kmid <- max(1, which(yu == mediany) - 1L)
  }
  # For large n, as.factor is slow
  # if(!is.factor(y)) y <- as.factor(y)
  if(is.factor(Y)) {
    ylev <- levels(Y)
    Y    <- unclass(Y)    # integer
  }
  else {
    ylev <- sort(unique(Y))
    Y    <- match(Y, ylev)
  }
  if(! ynumeric) {
    mediany <- quantile(Y, probs=.5, type=1L)
    # Find intercept close to median Y
    kmid    <- max(1, which(1L : length(ylev) == mediany) - 1L)
  }

  numy        <- tabulate(Y)
  names(numy) <- ylev

  k <- length(ylev)
  
  pppo <- if(length(ppo)) ncol(Z) else 0

  nrp <- length(ylev) - 1L
  if(nrp == 1) kmid <- 1

	ass     <- DesignAssign(atr, nrp, Terms)
  
  priorsd <- rep(priorsd, length=p)
	d <- list(X=if(! length(time)) Xs,
            y=if(! length(time)) Y,
            N=n, p=p, k=k, conc=conc,
            sds=as.array(priorsd))

  Nc <- 0
  if(length(cluster)) {
    d$psigma  <- psigma
    d$rsdmean <- rsdmean
    d$rsdsd   <- rsdsd
    cl        <- as.integer(as.factor(cluster))
    Nc        <- max(cl, na.rm=TRUE)
    d$Nc      <- Nc
    if(! length(time)) d$cluster <- cl  #for AR(1) cluster=rows of X,y
    method    <- 'sampling'
  }

  if(length(time)) {
    if(! all(time == floor(time))) stop('aTime variable must be integer-valued')
    d$N       <- NULL
    tim       <- time - min(time) + 1    # start at 1
    Nt        <- max(tim)
    Ntobs     <- length(unique(tim))
    d$Nt      <- Nt

    ## Reduce X matrix to first row (min time) per cluster
    first  <- tapply(1 : n, cl, function(i) i[which.min(tim[i])])
    Xbase  <- X[first,, drop=FALSE]
    wqrX   <- selectedQr(Xbase, center=TRUE, not=notransX)
    Xs     <- wqrX$X
    xbar   <- wqrX$xbar
    xqrsd <- apply(Xs, 2, sd)
    wqrX$X <- NULL

    ## Create integer ordinal response a Nc x Nt matrix
    ## Initially population with zeros, which will remain for
    ## missing assessments
    Yint <- matrix(0L, nrow=Nc, ncol=Nt)
    Yint[cbind(cl, tim)] <- Y
    d$X <- Xs
    d$y <- Yint
  }

  if(length(ppo)) {
    d$Z <- Zs
    d$q <- ncol(Z)
    priorsdppo <- rep(priorsdppo, length=pppo)
    d$sdsppo   <- as.array(priorsdppo)
  }
  
  if(standata) return(d)

	if(any(is.na(Xs)) | any(is.na(Y))) stop('program logic error')
  stanloc <- .Options$stancompiled
  if(! length(stanloc)) stop('options(stancompiled) not defined')
  
  fitter <- if(length(ppo)) ifelse(length(cluster), 'lrmcppo', 'lrmppo')
            else
              if(length(cluster) == 0) 'lrm'
            else
              if(length(time)) 'lrmcar1'
            else
              'lrmc'
  
  mfile <- paste0(stanloc, '/', fitter, '.rds')
  if(! file.exists(mfile))
    stop('you did not run rms::stanCompile to compile Stan code')
  mod      <- readRDS(mfile)
  stancode <- rstan::get_stancode(mod)

  ## See if previous fit had identical inputs and Stan code
  hashobj  <- list(d, inito, inits, iter, chains, loo, ppairs, method,
                   stancode, ...)
  datahash <- digest::digest(hashobj)
  if(length(prevhash) && prevhash == datahash) return(prevfit)
  
  hashobj  <- NULL
  
  itfailed <- function(w) is.list(w) && length(w$fail) && w$fail

  opt <- parm <- taus <- NULL

  if(method != 'sampling') {
    # Temporarily make concentration parameter = 1.0
    d$conc <- 1.0
    
    otime <- system.time(g <- rstan::optimizing(mod, data=d, init=inito))
    d$conc <- conc   # restore
    if(g$return_code != 0) {
      warning(paste('Optimizing did not work; return code', g$return_code,
                    '\nPosterior modes not computed'))
      opt <- list(coefficients=NULL,
                  sigmag=NA, deviance=NA,
                  return_code=g$return_code, hessian=NULL,
                  executionTime=otime)
    } else {
      parm <- g$par
      nam <- names(parm)
      al  <- nam[grep('alpha\\[', nam)]
      be  <- nam[grep('beta\\[',  nam)]
      ta  <- nam[grep('tau\\[',   nam)]
      alphas <- parm[al]
      betas  <- parm[be]
      betas <- matrix(betas, nrow=1) %*% t(wqrX$Rinv)
      names(betas)  <- atr$colnames
      alphas <- alphas - sum(betas * xbar)
      if(length(ppo)) {
        taus  <- matrix(parm[ta], nrow=pppo, ncol=k-2)
        taus  <- wqrZ$Rinv %*% taus
        alphas[-1] <- alphas[-1] - matrix(zbar, ncol=pppo) %*% taus
        ro     <- as.integer(gsub('tau\\[(.*),.*',    '\\1', ta))  # y cutoff
        co     <- as.integer(gsub('tau\\[.*,(.*)\\]', '\\1', ta))  # Z column
        namtau <- paste0(colnames(Z)[ro], ':y>=', ylev[-(1:2)][co])
        names(taus) <- namtau
      }
      names(alphas) <- if(nrp == 1) 'Intercept' else paste0('y>=', ylev[-1])
      
      opt <- list(coefficients=c(alphas, betas, taus),
                  sigmag=parm['sigmag'], deviance=-2 * g$value,
                  return_code=g$return_code, hessian=g$hessian,
                  executionTime=otime)
      }
    if(method == 'optimizing') return(opt)
  }

  if(progress != '') sink(progress, append=TRUE)
  
  #init <- NULL
  # if(length(parm)) {
  #   nam <- names(parm)
  #   nam <- nam[c(grep('alpha', nam), grep('beta', nam), grep('omega', nam),
  #                grep('pi', nam), grep('tau', nam))]  #, grep('theta', nam))]
  #   parm <- as.list(parm[nam])
  #   init <- function() parm
  #}

  debug(1)
  exclude <- c('sigmaw', 'gamma_raw', 'pi')  #, 'OR')
  if(! loo) exclude <- c(exclude, 'log_lik')
  g <- rstan::sampling(mod, pars=exclude, include=FALSE,
                       data=d, iter=iter, chains=chains, refresh=refresh,
                       init=inits, ...)
  debug(2)
  
  if(progress != '') sink()
	nam <- names(g)
	al  <- nam[grep('alpha\\[', nam)]
	be  <- nam[grep('beta\\[', nam)]
  ga  <- nam[grep('gamma\\[', nam)]

  draws  <- as.matrix(g)
  ndraws <- nrow(draws)
	alphas <- draws[, al, drop=FALSE]
	betas  <- draws[, be, drop=FALSE]
  betas  <- betas %*% t(wqrX$Rinv)

  omega   <- NULL      # non-intercepts, non-slopes
  clparm  <- character(0)
  gammas  <- NULL
  
  if(length(cluster)) {
    omega   <- cbind(sigmag = draws[, 'sigmag'])
    clparm  <- 'sigmag'
    cle     <- draws[, ga, drop=FALSE]
    gammas  <- apply(cle, 2, median)      # posterior median per subject
  }

  debug(3)
  tau <- ta <- tauInfo <- NULL
  if(length(ppo)) {
    ta     <- nam[grep('tau\\[', nam)]
    clparm <- c(clparm, ta)
    ro     <- as.integer(gsub('tau\\[(.*),.*',    '\\1', ta))  # y cutoff
    co     <- as.integer(gsub('tau\\[.*,(.*)\\]', '\\1', ta))  # Z column
    xt     <- colnames(Z)[ro]
    yt     <- ylev[-(1:2)][co]
    namtau <- paste0(xt, ':y>=', yt)
    taus   <- draws[, ta, drop=FALSE]
    tauInfo <- data.frame(intercept=1 + co, name=namtau, x=xt, y=yt)
    ## Need taus as a 3-dim array for intercept correction for centering
    ## Also helps with undoing the QR orthogonalization
    mtaus      <- array(taus, dim=c(ndraws, pppo, k - 2))

    for(j in 1 : (k - 2))
      mtaus[, , j] <- dropadim(mtaus, 3, j) %*% t(wqrZ$Rinv)
    taus <- matrix(as.vector(mtaus), nrow=ndraws, ncol=pppo * (k -2))
    colnames(taus) <- namtau

#    zalphacorr <- sweep(mtaus, 2, zbar, '*')
#    dim(zalphacorr) <- dim(zalphacorr)[-2]  # collapse to matrix
#    zalphacorr <- matrix(zbar, ncol=pppo) %*% mtaus
#          alphas[-1] <- alphas[-1] - matrix(zbar, ncol=pppo) %*% taus
#    zalphacorr <- cbind(0, rowSums(sweep(mtaus, 2, zbar, '*'), dims=2))
    }

  debug(4)
  
  epsmed <- NULL
  if(length(time)) {
    omega   <- cbind(omega, rho=draws[, 'rho'])
    clparm  <- c(clparm, 'rho')
    ep      <- nam[grep('eps', nam)]
    ## Summarize AR(1) random effects by taking median over draws
    ro <- as.integer(gsub('eps_raw\\[(.*),.*',    '\\1', ep))  # cluster
    co <- as.integer(gsub('eps_raw\\[.*,(.*)\\]', '\\1', ep))  # time
    ## For each cluster and time compute median over draws
    uro <- unique(ro)
    uco <- unique(co)
    eps <- draws[, ep, drop=FALSE]
    epsmed <- matrix(NA, nrow=length(uro), ncol=length(uco))
    for(i in uro) for(j in uco)
      epsmed[i, j] <- median(eps[, ro==i & co==j])
    }

  debug(5)
  diagnostics <-
		tryCatch(rstan::summary(g, pars=c(al, be, clparm), probs=NULL),
             error=function(...) list(fail=TRUE))
  if(itfailed(diagnostics)) {
    warning('rstan::summary failed; see fit component diagnostics')
    diagnostics <- list(pars=c(al, be, clparm), failed=TRUE)
  } else diagnostics <- diagnostics$summary[,c('n_eff', 'Rhat')]
  
  debug(6)
  if(length(ppairs)) {
    if(is.character(ppairs)) png(ppairs, width=1000, height=1000, pointsize=11)
    pairs(g, pars=c(al, be, clparm))
    if(is.character(ppairs)) dev.off()
    }

  debug(7)
	# Back-scale to original data scale
	alphacorr <- rowSums(sweep(betas, 2, xbar, '*')) # was xbar/xsd
	alphas    <- sweep(alphas, 1, alphacorr, '-')
  if(length(ppo))
    for(i in 1 : ndraws)
      for(j in 3 : k)
        alphas[i, j - 1] <- alphas[i, j - 1, drop=FALSE] -
                            sum(mtaus[i, , j-2] * zbar)
                    
  #  alphas[, -1] <- alphas[, -1, drop=FALSE] - zalphacorr
	colnames(alphas) <- if(nrp == 1) 'Intercept' else paste0('y>=', ylev[-1])
	colnames(betas)  <- atr$colnames

	draws            <- cbind(alphas, betas, taus)

  debug(8)
  param <- rbind(mean=colMeans(draws), median=apply(draws, 2, median))
  if(method != 'sampling') {
    param <- rbind(mode=opt$coefficients, param)
    opt$coefficients <- NULL
  }

  Loo <- lootime <- NULL
  if(loo) {
    lootime <- system.time(
      Loo <- tryCatch(rstan::loo(g), error=function(...) list(fail=TRUE)))
    if(itfailed(Loo)) {
      warning('loo failed; try running on loo(stanGet(fit object)) for more information')
      Loo <- NULL
      }
  }

  debug(9)
  
	res <- list(call=call, fitter=fitter, link='logistic',
							draws=draws, omega=omega,
              gammas=gammas, eps=epsmed,
              param=param, priorsd=priorsd, priorsdppo=priorsdppo,
              psigma=psigma, rsdmean=rsdmean, rsdsd=rsdsd, conc=conc,
              notransX=notransXn, notransZ=notransZn, xqrsd=xqrsd,
              N=n, p=p, pppo=pppo, yname=yname, yunique=ylev, freq=numy,
						  alphas=al, betas=be, taus=ta, tauInfo=tauInfo,
						  xbar=xbar, Design=atr, scale.pred=c('log odds', 'Odds Ratio'),
							terms=Terms, assign=ass, na.action=atrx$na.action, fail=FALSE,
							non.slopes=nrp, interceptRef=kmid, sformula=sformula,
							x=if(x) X, y=if(y) Y, z=if(x) Z, loo=Loo, lootime=lootime,
              clusterInfo=if(length(cluster))
                list(cluster=if(x) cluster else NULL, n=Nc, name=clustername),
              timeInfo=if(length(time))
                         list(time=if(x) time else NULL, n=Nt, nobs=Ntobs,
                              name=timename),
							opt=opt, diagnostics=diagnostics,
              iter=iter, chains=chains, stancode=stancode, datahash=datahash)
	class(res) <- c('blrm', 'rmsb', 'rms')
  if(length(file)) saveRDS(res, file)
  res$rstan <- g
	res
}

##' Compute Indexes of Predictive Accuracy and Their Uncertainties
##'
##' For a binary or ordinal logistic regression fit from \code{blrm}, computes several indexes of predictive accuracy along with highest posterior density intervals for them.  Optionally plots their posterior densities.
##' When there are more than two levels of the outcome variable, computes Somers' Dxy and c-index on a random sample of 10,000 observations.
##' @title blrmStats
##' @param fit an object produced by \code{blrm}
##' @param ns number of posterior draws to use in the calculations (default is 400)
##' @param prob HPD interval probability (default is 0.95)
##' @param pl set to \code{TRUE} to plot the posterior densities using base graphics
##' @param dist if \code{pl} is \code{TRUE} specifies whether to plot the density estimate (the default) or a histogram
##' @return list of class \code{'blrmStats'} whose most important element is \code{Stats}.  The indexes computed are defined below, with gp, B, EV, and vp computed using the intercept corresponding to the median value of Y.  See \url{https://fharrell.com/post/addvalue} for more information.
##' \describe{
##'  \item{"Dxy"}{Somers' Dxy rank correlation between predicted and observed.  The concordance probability (c-index; AUROC in the binary Y case) may be obtained from the relationship Dxy=2(c-0.5).}
##'  \item{"g"}{Gini's mean difference: the average absolute difference over all pairs of linear predictor values}
##'  \item{"gp"}{Gini's mean difference on the predicted probability scale}
##'  \item{"B"}{Brier score}
##'  \item{"EV"}{explained variation}
##'  \item{"v"}{variance of linear predictor}
##'  \item{"vp"}{variable of estimated probabilities}
##' }
##' @seealso [Hmisc::rcorr.cens]
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   blrmStats(f, pl=TRUE)   # print and plot
##' }
##' @author Frank Harrell
blrmStats <- function(fit, ns=400, prob=0.95, pl=FALSE,
                      dist=c('density', 'hist')) {
  dist <- match.arg(dist)
  
  f <- fit[c('x', 'y', 'z', 'non.slopes', 'interceptRef', 'pppo',
             'draws', 'yunique', 'tauInfo')]
  X <- f$x
  Z <- f$z
  y <- f$y
  if(length(X) == 0 | length(y) == 0)
    stop('must have specified x=TRUE, y=TRUE to blrm')
  y <- as.integer(y) - 1
  nrp    <- f$non.slopes
  kmid   <- f$interceptRef  # intercept close to median
  s      <- tauFetch(f, intercept=kmid, what='nontau')
  # old fit objects did not include pppo for non-PPO models:
  pppo   <- if(! length(f$pppo)) 0 else f$pppo  
  if(pppo > 0) stau <-tauFetch(f, intercept=kmid, what='tau')
  ndraws <- nrow(s)
  ns     <- min(ndraws, ns)
  if(ns < ndraws) {
    j <- sample(1 : ndraws, ns, replace=FALSE)
    s <- s[j,, drop=FALSE]
    if(pppo > 0) stau <- stau[j,, drop=FALSE]
    }
  ylev  <- f$yunique
  ybin  <- length(ylev) == 2
  stats <- matrix(NA, nrow=ns, ncol=8)
  colnames(stats) <- c('Dxy', 'C', 'g', 'gp', 'B', 'EV', 'v', 'vp')
  dxy <- if(length(ylev) == 2)
           function(x, y) somers2(x, y)['Dxy']
         else
           function(x, y) {
             con <- survival::survConcordance.fit(Surv(y), x)
             conc <- con['concordant']; disc <- con['discordant']
             - (conc - disc) / (conc + disc)
             }
  brier <- function(x, y) mean((x - y) ^ 2)
  br2   <- function(p) var(p) / (var(p) + sum(p * (1 - p)) / length(p))

  nobs <- length(y)
  is <- if((nobs <= 10000) || (length(ylev) == 2)) 1 : nobs else
              sample(1 : nobs, 10000, replace=FALSE)

  for(i in 1 : ns) {
    beta  <- s[i,, drop=FALSE ]
    lp    <- cbind(1, X) %*% t(beta)
    if(pppo > 0) {
      tau <- stau[i,, drop=FALSE]
      lp  <- lp + Z %*% t(tau)
      }
    prb  <- plogis(lp)
    d    <- dxy(lp[is], y[is])
    C    <- (d + 1.) / 2.
    st <- c(d, C, GiniMd(lp), GiniMd(prb),
            brier(prb, y > kmid - 1), 
            br2(prb), var(lp), var(prb))
    stats[i,] <- st
  }
  sbar  <- colMeans(stats)
  se    <- apply(stats, 2, sd)
  hpd   <- apply(stats, 2, HPDint, prob=prob)
  sym   <- apply(stats, 2, distSym)
  
  text <- paste0(round(sbar, 3), ' [', round(hpd[1,], 3), ', ',
                 round(hpd[2,], 3), ']')
  Stats <- rbind(Mean=sbar, SE=se, hpd=hpd, Symmetry=sym)
  statnames <- colnames(Stats)
  names(text) <- statnames

  if(pl) {
    par(mfrow=c(4,2), mar=c(3, 2, 0.5, 0.5), mgp=c(1.75, .55, 0))
    for(w in setdiff(statnames, 'C')) {
      p <- switch(dist,
             density = {
               den <- density(stats[, w])
               plot(den, xlab=w, ylab='', type='l', main='') },
             hist = {
               den <- hist(stats[, w], probability=TRUE,
                           nclass=100, xlab=w, ylab='', main='')
               den$x <- den$breaks; den$y <- den$density } )
      ref <- c(sbar[w], hpd[1, w], hpd[2, w])
      abline(v=ref, col=gray(0.85))
      text(min(den$x), max(den$y), paste0('Symmetry:', round(sym[w], 2)),
           adj=c(0, 1), cex=0.7)
    }
  }
  structure(list(stats=Stats, text=text, ndraws=ndraws, ns=ns,
                 non.slopes=nrp, intercept=kmid), class='blrmStats')
  }

##' Print Details for \code{blrmStats} Predictive Accuracy Measures
##'
##' Prints results of \code{blrmStats} with brief explanations
##' @title print.blrmStats
##' @param x an object produced by \code{blrmStats}
##' @param dec number of digits to round indexes
##' @param ... ignored
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   s <- blrmStats(...)
##'   s    # print with defaults
##'   print(s, dec=4)
##' }
##' @author Frank Harrell
print.blrmStats <- function(x, dec=3, ...) {
  ns <- x$ns; ndraws <- x$ndraws
  if(ns < ndraws)
    cat('Indexes computed for a random sample of', ns, 'of', ndraws,
        'posterior draws\n\n')
  else
    cat('Indexes computed on', ndraws, 'posterior draws\n\n')
  
  if(x$non.slopes > 1)
    cat('gp, B, EV, and vp are for intercept', x$intercept,
        'out of', x$non.slopes, 'intercepts\n\n')
  print(round(x$stats, dec))
  cat('\nDxy: 2*(C - 0.5)   C: concordance probability',
      'g: Gini mean |difference| on linear predictor (lp)',
      'gp: Gini on predicted probability        B: Brier score',
      'EV: explained variation on prob. scale   v: var(lp)   vp: var(prob)\n',
      sep='\n')
    }


##' Print \code{blrm} Results
##'
##' Prints main results from \code{blrm} along with indexes and predictive accuracy and their highest posterior density intervals computed from \code{blrmStats}.
##' @title print.blrm
##' @param x object created by \code{blrm}
##' @param dec number of digits to print to the right of the decimal
##' @param coefs specify \code{FALSE} to suppress printing parameter estimates, and in integer k to print only the first k
##' @param intercepts set to \code{FALSE} to suppress printing intercepts.   Default is to print them unless there are more than 9.
##' @param prob HPD interval probability for summary indexes
##' @param ns number of random samples of the posterior draws for use in computing HPD intervals for accuracy indexes
##' @param title title of output
##' @param ... passed to \code{prModFit}
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   options(lang='html')   # default is lang='plain'; also can be latex
##'   f               # print using defaults
##'   print(f, posterior.summary='median')   # instead of post. means
##' }
##' @author Frank Harrell
print.blrm <- function(x, dec=4, coefs=TRUE, intercepts=x$non.slopes < 10,
                       prob=0.95, ns=400,
                       title='Bayesian Logistic Regression Model', ...) {
  latex <- prType() == 'latex'
  
  z <- list()
  k <- 0
  
  if(length(x$freq) > 3 && length(x$freq) < 50) {
    k <- k + 1
    z[[k]] <- list(type='print', list(x$freq),
                   title='Frequencies of Responses')
  }
  if(length(x$na.action)) {
    k <- k + 1
    z[[k]] <- list(type=paste('naprint',class(x$na.action),sep='.'),
                   list(x$na.action))
  }

  qro <- function(x) {
    r <- round(c(median(x), HPDint(x, prob)), 4)
    paste0(r[1], ' [', r[2], ', ', r[3], ']')
    }
  ci  <- x$clusterInfo
  sigmasum <- NULL
  if(length(ci)) sigmasum <- qro(x$omega[, 'sigmag'])

  ti     <- x$timeInfo
  rhosum <- if(length(ti)) qro(x$omega[, 'rho'])
 
  loo <- x$loo
  elpd_loo <- p_loo <- looic <- NULL
  if(length(loo)) {
    lo <- loo$estimates
    pm <- if(prType() == 'plain') '+/-' else
              markupSpecs[[prType()]][['plminus']]
    nlo <- rownames(lo)
     lo <- paste0(round(lo[, 'Estimate'], 2), pm, round(lo[, 'SE'], 2))
    elpd_loo <- lo[1]; p_loo <- lo[2]; looic <- lo[3]
    }
  misc <- reListclean(Obs             = x$N,
                      Draws           = nrow(x$draws),
                      Chains          = x$chains,
                      Imputations     = x$n.impute,
                      p               = x$p,
                      'Cluster on'    = ci$name,
                      Clusters        = ci$n,
                      'sigma gamma'   = sigmasum,
                      'Time variable' = ti$name,
                      'Possible times'= ti$n,
                      'Observed times'= ti$nobs,
                      rho             = rhosum)
  
  if(length(x$freq) < 4) {
    names(x$freq) <- paste(if(latex)'~~' else ' ',
                           names(x$freq), sep='')
    misc <- c(misc[1], x$freq, misc[-1])
  }
  
  a <- blrmStats(x, ns=ns)$text

  mixed <- reListclean('LOO log L'  = elpd_loo,
                      'LOO IC'      = looic,
                      'Effective p' = p_loo,
                      B             = a['B'])
  
  disc <- reListclean(g       = a['g'],
                      gp      = a['gp'],
                      EV      = a['EV'],
                      v       = a['v'],
                      vp      = a['vp'])
                     
  discr <-reListclean(C       = a['C'],
                      Dxy     = a['Dxy'])

  
  headings <- c('','Mixed Calibration/\nDiscrimination Indexes',
                   'Discrimination\nIndexes',
                   'Rank Discrim.\nIndexes')
  
  data <- list(misc, c(mixed, NA), c(disc, NA), c(discr, NA))
  k <- k + 1
  z[[k]] <- list(type='stats', list(headings=headings, data=data))

  if(coefs) {
    k <- k + 1
    z[[k]] <- list(type='coefmatrix',
                   list(bayes=print.rmsb(x, prob=prob, intercepts=intercepts,
                                         pr=FALSE)))
  }

  footer <- if(length(x$notransX))
              paste('The following parameters remained separate (where not orthogonalized) during model fitting so that prior distributions could be focused explicitly on them:',
                    paste(x$notransX, collapse=', '))

  prModFit(x, title=title, z, digits=dec, coefs=coefs, footer=footer, ...)
}

##' Make predictions from a \code{blrm} fit
##'
##' Predict method for \code{blrm} objects
##' @title predict.blrm
##' @param object,...,type,se.fit,codes see [predict.lrm] 
##' @param posterior.summary set to \code{'median'} or \code{'mode'} to use posterior median/mode instead of mean
##' @param cint probability for highest posterior density interval
##' @return a data frame,  matrix, or vector with posterior summaries for the requested quantity, plus an attribute \code{'draws'} that has all the posterior draws for that quantity.  For \code{type='fitted'} and \code{type='fitted.ind'} this attribute is a 3-dimensional array representing draws x observations generating predictions x levels of Y.
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   predict(f, newdata, type='...', posterior.summary='median')
##' }
##' @seealso [predict.lrm]
##' @author Frank Harrell
##' @md
predict.blrm <- function(object, ..., 
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		  "terms", "cterms", "ccterms", "adjto", "adjto.data.frame",
      "model.frame"),
		se.fit=FALSE, codes=FALSE,
    posterior.summary=c('mean', 'median'), cint=0.95) {
  
  type              <- match.arg(type)
  posterior.summary <- match.arg(posterior.summary)

  if(se.fit) warning('se.fit does not apply to Bayesian models')

  if(object$pppo > 0 && (type %in% c('lp', 'mean', 'fitted', 'fitted.ind')))
    stop('type not yet implemented for partial proportional odds models')

  if(type %in% c('lp', 'x', 'data.frame', 'terms', 'cterms', 'ccterms',
                 'adjto', 'adjto.data.frame', 'model.frame'))
    return(predictrms(object,...,type=type,
                      posterior.summary=posterior.summary))

  if(type %nin% c('mean', 'fitted', 'fitted.ind'))
    stop('types other than mean, fitted, fitted.ind not yet implemented')

  ns      <- object$non.slopes
  tauinfo <- object$tauInfo

  if(ns == 1 & type == "mean")
    stop('type="mean" makes no sense with a binary response')

  draws   <- object$draws
  ndraws  <- nrow(draws)
  cn      <- colnames(draws)
  ints    <- draws[, 1 : ns, drop=FALSE]
  betanam <- setdiff(cn, c(cn[1:ns], tauinfo$name))
  betas   <- draws[, betanam, drop=FALSE]
  
  X <- predictrms(object, ..., type='x', posterior.summary=posterior.summary)
  rnam  <- rownames(X)
  n     <- nrow(X)

  ## Get cumulative probability function used
  link  <- object$link
  if(! length(link)) link <- 'logistic'
  cumprob <- probabilityFamilies[[link]]$cumprob
  
  if(ns == 1) return(cumprob(ints + betas %*% t(X)))

  cnam  <- cn[1:ns]
  ylev  <- object$yunique
  if(type == 'mean') {
    if(codes) vals <- 1:length(ylev)
    else {
      vals <- as.numeric(ylev)
      if(any(is.na(vals)))
         stop('values of response levels must be numeric for type="mean" and codes=FALSE')
    }
  }

  ynam <- paste(object$yname, "=", ylev, sep="")
  PP   <- array(NA, dim=c(ndraws, n, ns),
                dimnames=list(NULL, rnam, cnam))
  PPeq <- array(NA, dim=c(ndraws, n, ns + 1),
                dimnames=list(NULL, rnam, ynam))

  M  <- matrix(NA, nrow=ndraws, ncol=n, dimnames=list(NULL, rnam))

  for(i in 1 : ndraws) {
    inti    <- ints[i,]
    betai   <- betas[i,,drop=FALSE]
    xb      <- betai %*% t(X)
    xb      <- sapply(inti, "+", xb)  # n x ns matrix, col j adds intercept j
    P       <- cumprob(xb)            # preserves matrix;  1 x nrow(X)
    PP[i,,] <- P
    if(type != 'fitted') {
      Peq       <- cbind(1, P) - cbind(P, 0)
      PPeq[i,,] <- Peq
      if(type == 'mean') M[i, ] <- drop(Peq %*% vals)
      }
  }

  ps <- switch(posterior.summary, mean=mean, median=median)
  h <- function(x) {
    s  <- ps(x)
    ci <- HPDint(x, cint)
    r  <- c(s, ci)
    names(r)[1] <- posterior.summary
    r
    }

  ## Function to summarize a 3-d array and transform the results
  ## to a data frame with variables x (row names of predictions), y level,
  ## mean/median, Lower, Upper
  ## Input is ndraws x # predicton requests x # y categories (less 1 if cum.p.)
  
  s3 <- function(x) {
    yl <- if(type == 'fitted') cnam else ynam
    d <- expand.grid(y = yl, x = rnam, stat=NA, Lower=NA, Upper=NA,
                     stringsAsFactors=FALSE)
    for(i in 1 : nrow(d)) {
      u <- h(x[, d$x[i], d$y[i]])
      d$stat[i] <- u[1]
      d$Lower[i] <- u['Lower']
      d$Upper[i] <- u['Upper']
    }
    names(d)[3] <- upFirst(posterior.summary)
    d
  }

  ## Similar for a 2-d array
  s2 <- function(x) {
    d <- expand.grid(x = rnam, stat=NA, Lower=NA, Upper=NA,
                     stringsAsFactors=FALSE)
    for(i in 1 : nrow(d)) {
      u <- h(x[, d$x[i]])
      d$stat[i]  <- u[1]
      d$Lower[i] <- u['Lower']
      d$Upper[i] <- u['Upper']
    }
    names(d)[2] <- upFirst(posterior.summary)
    d
  }
  
  
  r <- switch(type,
              fitted     =  structure(s3(PP), draws=PP),
              fitted.ind =  structure(s3(PPeq), draws=PPeq),
              mean       =  structure(s2(M), draws=M))

  class(r) <- c('predict.blrm', class(r))
  r
}

##' Print Predictions for \code{blrm}
##'
##' Prints the summary portion of the results of \code{predict.blrm}
##' @title print.predict.blrm
##' @param x result from \code{predict.blrm}
##' @param digits number of digits to round numeric results
##' @param ... ignored
##' @author Frank Harrell
print.predict.blrm <- function(x, digits=3, ...) {
  numvar <- sapply(x, is.numeric)
  numvar <- names(numvar)[numvar]
  for(j in numvar) x[, j] <- round(x[, j], digits)
  print.data.frame(x)
}


##' Function Generator for Mean Y for \code{blrm}
##'
##' Creates a function to turn a posterior summarized linear predictor lp (e.g. using posterior mean of intercepts and slopes) computed at the reference intercept into e.g. an estimate of mean Y using the posterior mean of all the intercepts
##' @title Mean.blrm
##' @param object a \code{blrm} fit
##' @param codes if \code{TRUE}, use the integer codes \eqn{1,2,\ldots,k} for the \eqn{k}-level response in computing the predicted mean response.
##' @param posterior.summary defaults to posterior mean; may also specify \code{"median"}.  Must be consistent with the summary used when creating \code{lp}.
##' @param ... ignored
##' @return an R function
##' @author Frank Harrell
Mean.blrm <- function(object, codes=FALSE,
                      posterior.summary=c('mean', 'median'), ...) {
  posterior.summary <- match.arg(posterior.summary)
  ns <- num.intercepts(object)
  if(ns < 2)
stop('using this function only makes sense for >2 ordered response categories')

  ## Get cumulative probability function used
  link  <- object$link
  if(! length(link)) link <- 'logistic'
  cumprob <- probabilityFamilies[[link]]$cumprob

  ylev   <- object$yunique
  if(codes) vals <- 1:length(ylev)
    else {
      vals <- as.numeric(ylev)
      if(any(is.na(vals)))
         stop('values of response levels must be numeric for type="mean" and codes=FALSE')
    }

  cof <- coef(object, stat=posterior.summary)
  intercepts <- cof[1 : ns]
    
  f <- function(lp=numeric(0), intercepts=numeric(0), values=numeric(0),
                  interceptRef=integer(0), cumprob=cumprob) {
    ns <- length(intercepts)
    lp <- lp - intercepts[interceptRef]
    xb <- sapply(intercepts, '+', lp)
    P  <- matrix(cumprob(xb), ncol=ns)
    P  <- cbind(1, P) - cbind(P, 0)
    m  <- drop(P %*% values)
    names(m) <- names(lp)
    m
  }
 
  ir <- object$interceptRef
  if(! length(ir)) ir <- 1
  formals(f) <- list(lp=numeric(0), intercepts=intercepts,
                     values=vals, interceptRef=ir, cumprob=cumprob)
  f
}

##' Function Generator for Quantiles of Y for \code{blrm}
##'
##' Creates a function to turn a posterior summarized linear predictor lp (e.g. using posterior mean of intercepts and slopes) computed at the reference intercept into e.g. an estimate of a quantile of Y using the posterior mean of all the intercepts
##' @title Quantile.blrm
##' @param object a \code{blrm} fit
##' @param codes if \code{TRUE}, use the integer codes \eqn{1,2,\ldots,k} for the \eqn{k}-level response in computing the quantile
##' @param posterior.summary defaults to posterior mean; may also specify \code{"median"}.  Must be consistent with the summary used when creating \code{lp}.
##' @param ... ignored
##' @return an R function
##' @author Frank Harrell
Quantile.blrm <- function(object, codes=FALSE,
                      posterior.summary=c('mean', 'median'), ...) {
  posterior.summary <- match.arg(posterior.summary)
  ns <- num.intercepts(object)
  if(ns < 2)
stop('using this function only makes sense for >2 ordered response categories')

  ## Get cumulative probability function used
  link  <- object$link
  if(! length(link)) link <- 'logistic'
  pfam <- probabilityFamilies[[link]]
  cumprob <- pfam$cumprob
  inverse <- pfam$inverse

  ylev   <- object$yunique
  if(codes) vals <- 1:length(ylev)
    else {
      vals <- as.numeric(ylev)
      if(any(is.na(vals)))
         stop('values of response levels must be numeric for type="mean" and codes=FALSE')
    }

  cof        <- coef(object, stat=posterior.summary)
  intercepts <- cof[1 : ns]

  f <- function(q=.5, lp=numeric(0), intercepts=numeric(0), values=numeric(0),
                  interceptRef=integer(0), cumprob=NULL, inverse=NULL) {
        ## Solve inverse(1 - q) = a + xb; inverse(1 - q) - xb = a
        ## Shift intercepts to left one position because quantile
        ## is such that Prob(Y <= y) = q whereas model is stated in
        ## Prob(Y >= y)
        lp <- lp - intercepts[interceptRef]
        ## Interpolation on linear predictors scale:
        ## z <- approx(c(intercepts, -1e100), values,
        ##             xout=inverse(1 - q) - lp, rule=2)$y
        ## Interpolation approximately on Y scale:
        lpm <- mean(lp, na.rm=TRUE)
        z <- approx(c(cumprob(intercepts + lpm), 0), values,
                    xout=cumprob(inverse(1 - q) - lp + lpm), rule=2)$y
        names(z) <- names(lp)
        z
      }
    trans <- object$trans
    formals(f) <- list(q=.5, lp=numeric(0), intercepts=intercepts,
                       values=vals, interceptRef=object$interceptRef,
                       cumprob=cumprob, inverse=inverse)
    f
}

##' Function Generator for Exceedance Probabilities for \code{blrm}
##'
##' For a \code{blrm} object generates a function for computing the	estimates of the function Prob(Y>=y) given one or more values of thelinear predictor using the reference (median) intercept.  This function can optionally be evaluated at only a set of user-specified	\code{y} values, otherwise a right-step function is returned.  There is a plot method for plotting the step functions, and if more than one linear predictor was evaluated multiple step functions are drawn. \code{ExProb} is especially useful for \code{\link{nomogram}}.  The linear predictor argument is a posterior summarized linear predictor lp (e.g. using posterior mean of intercepts and slopes) computed at the reference intercept
##' @title ExProb.blrm
##' @param object a \code{blrm} fit
##' @param codes if \code{TRUE}, use the integer codes \eqn{1,2,\ldots,k} for the \eqn{k}-level response in computing the quantile
##' @param posterior.summary defaults to posterior mean; may also specify \code{"median"}.  Must be consistent with the summary used when creating \code{lp}.
##' @param ... ignored
##' @return an R function
##' @author Frank Harrell
ExProb.blrm <- function(object, codes=FALSE,
                        posterior.summary=c('mean', 'median'), ...) {
  posterior.summary <- match.arg(posterior.summary)
  ns <- num.intercepts(object)
  if(ns < 2)
stop('using this function only makes sense for >2 ordered response categories')

  ## Get cumulative probability function used
  link  <- object$link
  if(! length(link)) link <- 'logistic'
  pfam <- probabilityFamilies[[link]]
  cumprob <- pfam$cumprob
  
  ylev   <- object$yunique
  if(codes) vals <- 1:length(ylev)
    else {
      vals <- as.numeric(ylev)
      if(any(is.na(vals)))
         stop('values of response levels must be numeric for type="mean" and codes=FALSE')
    }

  cof        <- coef(object, stat=posterior.summary)
  intercepts <- cof[1 : ns]

  f <- function(lp=numeric(0), y=NULL, intercepts=numeric(0),
                values=numeric(0),
                interceptRef=integer(0), cumprob=NULL, yname=NULL) {
    lp <- lp - intercepts[interceptRef]
    prob <- cumprob(sapply(c(1e30, intercepts), '+', lp))
    dim(prob) <- c(length(lp), length(values))
        
    if(! length(y)) {
      colnames(prob) <- paste('Prob(Y>=', values, ')', sep='')
      return(structure(list(y=values, prob=prob, yname=yname),
                       class='ExProb'))
    }
    p <- t(apply(prob, 1,
                 function(probs) {
                   pa <- approx(values, probs, xout=y,
                                f=1, method='constant')$y
                   pa[y < min(values)] <- 1.
                   pa[y > max(values)] <- 0.
                   pa
                 } ) )
    if(length(y) > 1) {
      colnames(p) <- paste('Prob(Y>=', y, ')', sep='')
      p <- list(y=y, prob=p)
    }
    structure(drop(p), class='ExProb')
  }
  trans <- object$trans
  formals(f) <- list(lp=numeric(0), y=NULL, intercepts=intercepts,
                     values=vals, interceptRef=object$interceptRef,
                     cumprob=cumprob, yname=object$yname)
  f
}

##' Fetch Partial Proportional Odds Parameters
##'
##' Fetches matrix of posterior draws for partial proportional odds parameters (taus) for a given intercept.  Can also form a matrix containing both regular parameters and taus, or for just non-taus.
##' @title tauFetch
##' @param fit an object created by \code{blrm}
##' @param intercept integer specifying which intercept to fetch
##' @param what specifies the result to return
##' @return matrix with number of raws equal to the numnber of original draws
##' @author Frank Harrell
tauFetch <- function(fit, intercept, what=c('tau', 'nontau', 'both')) {
  what   <- match.arg(what)
  f      <- fit[c('tauInfo', 'draws', 'non.slopes')]
  info   <- f$tauInfo
  draws  <- f$draws
  nd     <- nrow(draws)
  cn     <- colnames(draws)
  int    <- intercept
  nints  <- f$non.slopes
  if(int > nints) stop('intercept is too large')
  ## Keep only the intercept of interest
  cn     <- c(cn[intercept], cn[-(1 : nints)])
  nontau <- setdiff(cn, info$name)
  i      <- if(int > 1) subset(info, intercept == int)$name
  # Partial proportional odds parameters start with the 2nd intercept
  switch(what,
         tau    = if(! length(i)) matrix(0, nrow=nd, ncol=1)
                  else draws[, i, drop=FALSE],
         nontau = draws[, nontau, drop=FALSE],
         both   = if(! length(i))
                    cbind(draws[, nontau, drop=FALSE], 0) else 
                              draws[, c(nontau, i), drop=FALSE]
        )
}
