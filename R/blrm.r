##' Bayesian Binary and Ordinal Logistic Regression
##'
##' Uses \code{rstan} with pre-compiled Stan code whose location is given by the user in \code{options(stancompiled='...')} to get posterior draws of parameters from a binary logistic or proportional odds semiparametric ordinal logistic model.  The Stan code internally using the qr decompositon on the design matrix so that highly collinear columns of the matrix do not hinder the posterior sampling.  The parameters are transformed back to the original scale before returning results to R.   Design matrix columns re centered before running Stan, so Stan diagnostic output will have the intercept terms shifted but the results of \code{blrm} for intercepts are for the original uncentered data.  The only prior distributions for regression betas are normal with mean zero, and the vector of prior standard deviations is given in \code{priorsd}.  These priors are for the qr-projected design matrix elements, except that the very last element is not changed.  So if one has a single non-interactive linear or binary variable for which a skeptical prior is designed, put that variable last in the model.
##'
##' \code{blrm} also handles single-level hierarchical random effects models for the case when there are repeated measurements per subject which are reflected as random intercepts, and a different model that allows for AR(1) serial correlation within subject.  For both setups, a \code{cluster} term in the model signals the existence of subject-specific random effects, and an additional model term \code{aTime(time variable)} signals the use of the AR(1) within-subject model.  The \code{time} variable must be integer valued and there can be arbitrary gaps between measurements.  However if the maximum time exceeds 200 or so one can expect much longer computation time.  When \code{aTime()} is present, the cluster-specific random effects then become the random effect for a subject's \code{time=1} record.  When \code{aTime} is specified, the covariates must not change over records within subject (no time-dependent covariates), as only the covariate vector for the earliest time for a subject is used.
##' @title blrm
##' @param formula a R formula object that can use \code{rms} package enhancements such as the restricted interaction operator
##' @param data a data frame
##' @param subset a logical vector or integer subscript vector specifying which subset of data whould be used
##' @param na.action default is \code{na.delete} to remove missings and report on them
##' @param priorsd vector of prior standard deviations.  If the vector is shorter than the number of model parameters, it will be repeated until the length equals the number of parametertimes.
##' @param rsdmean the assumed mean of the prior distribution of the standard deviation of random effects.  An exponential prior distribution is assumed, and the rate for that distribution is the reciprocal of the mean.  The default is a mean of 1.0, which is reasonable for a unitless regression model scale such as log odds.
##' @param ar1sdmean the assumed mean of the prior distribution of the standard deviation of within-subject white noise.   The setup is the same as with \code{rsdmean}.
##' @param iter number of posterior samples per chain for [rstan::sampling] to run
##' @param chains number of separate chains to run
##' @param refresh see [rstan::sampling]
##' @param x set to \code{FALSE} to not store the design matrix in the fit.  \code{x=TRUE} is needed if running \code{blrmStats} for example.
##' @param y set to \code{FALSE} to not store the response variable in the fit
##' @param loo set to \code{FALSE} to not run \code{loo} and store its result as object \code{loo} in the returned object
##' @param method set to \code{'optimizing'} to run the Stan optimizer and not do posterior sampling, \code{'both'} (the default) to run both the optimizer and posterior sampling, or \code{'sampling'} to run only the posterior sampling and not compute posterior modes. Running \code{optimizing} is a way to obtain maximum likelihood estimates and allows one to quickly study the effect of changing the prior distributions.  When \code{method='optimizing'} is used the result returned is not a standard \code{blrm} object but is instead the parameter estimates, -2 log likelihood, and optionally the Hession matrix (if you specify \code{hessian=TRUE} in ...).  When \code{method='both'} is used, \code{rstan::sampling} and \code{rstan::optimizing} are both run, and parameter estimates (posterior modes) from \code{optimizing} are stored in a matrix \code{param} in the fit object, which also contains the posterior means and medians, and other results from \code{optimizing} are stored in object \code{opt} in the \code{blrm} fit object.  When random effects are present, \code{method} is automatically set to \code{'sampling'} as maximum likelihood estimates without marginalizing over the random effects do not make sense.
##' @param standata set to \code{TRUE} to return the Stan data list and not run the model
##' @param ... passed to \code{rstan::sampling} or \code{rstan:optimizing}
##' @return an \code{rms} fit object of class \code{blrm}, \code{rmsb}, \code{rms} that also contains \code{rstan} results under the name \code{rstan}.  In the \code{rstan} results, which are also used to produce diagnostics, the intercepts are shifted because of the centering of columns of the design matrix done by \code{blrm}.  With \code{method='optimizing'} a class-less list is return with these elements: \code{coefficients} (MLEs), \code{theta} (non-intercept parameters on the QR decomposition scale), \code{deviance} (-2 log likelihood), \code{return_code} (see \code{rstan::optimizing}), and, if you specified \code{hessian=TRUE} to \code{blrm}, the Hessian matrix.
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
##'   traceplot(s)        # Stan diagnostic plots by chain
##'   traceplot(s, pars=f$betas)  # Same but only for beta parameters
##'   blrmStats(f)        # more details about predictive accuracy measures
##'   ggplot(Predict(...))   # standard rms output
##'   summary(f, ...)     # invokes summary.rms
##'   contrast(f, ...)    # contrast.rms computes credible intervals
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
blrm <- function(formula, data, subset, na.action=na.delete,
								 priorsd=rep(100, p), rsdmean=1, ar1sdmean=1,
								 iter=2000, chains=4, refresh=0,
								 x=TRUE, y=TRUE, loo=TRUE,
                 method=c('both', 'sampling', 'optimizing'),
                 standata=FALSE,
                 ...) {

	call <- match.call()
  m <- match.call(expand.dots=FALSE)
  mc <- match(c("formula", "data", "subset", "na.action"), 
             names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action
  m$drop.unused.levels <- TRUE

  m[[1]] <- as.name("model.frame")
  nact <- NULL
  if(missing(data)) data <- NULL

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

  method <- match.arg(method)

  X <- Design(eval.parent(m))   # Design handles cluster()
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

  if(! length(time)) {
    Xs  <- scale(X, center=TRUE, scale=FALSE)
    ## scinfo <- attributes(Xs)[c('scaled:center', 'scaled:scale')]
    ## xbar   <- as.vector(scinfo[[1]])
    ## xsd    <- as.vector(scinfo[[2]])
    xbar   <- as.vector(attr(Xs, 'scaled:center'))
    ## if(any(xsd == 0)) stop('a variable is constant')
    }
	
	n <- nrow(X)
	p <- ncol(X)
	Y <- as.factor(Y)
	ylev <- levels(Y)
	yint <- as.integer(Y)
	k    <- length(ylev)

  ## Find intercept that is close to the median of y
  mediany <- quantile(yint, probs=.5, type=1L)
  kmid    <- max(1, which(1L : length(ylev) == mediany) - 1L)

  nrp <- length(ylev) - 1L
  if(nrp == 1) kmid <- 1

	ass <- DesignAssign(atr, nrp, Terms)
  priorsd <- rep(priorsd, length=p)
	d <- list(X=if(! length(time)) Xs,
            y=if(! length(time)) yint,
            N=n, p=p, k=k, sds=as.array(priorsd),
            rate = 1. / rsdmean)
  Nc <- 0
  if(length(cluster)) {
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
    d$ratew   <- 1. / ar1sdmean

    ## Reduce X matrix to first row (min time) per cluster
    first <- tapply(1 : n, cl, function(i) i[which.min(tim[i])])
    Xbase <- X[first,, drop=FALSE]
    Xs    <- scale(Xbase, center=TRUE, scale=FALSE)
    xbar  <- as.vector(attr(Xs, 'scaled:center'))

    ## Create integer ordinal response a Nc x Nt matrix
    ## Initially population with zeros, which will remain for
    ## missing assessments
    Yint <- matrix(0L, nrow=Nc, ncol=Nt)
    Yint[cbind(cl, tim)] <- yint
    d$X <- Xs
    d$y <- Yint
    }
  if(standata) return(d)
  
	if(any(is.na(Xs)) | any(is.na(yint))) stop('program logic error')
  stanloc <- .Options$stancompiled
  if(! length(stanloc)) stop('options(stancompiled) not defined')
  
  fitter <- if(length(cluster) == 0) 'lrmqr'
            else
              if(length(time)) 'lrmqrcar1'
            else
              'lrmqrc'
  file <- paste0(stanloc, '/', fitter, '.rds')
  if(! file.exists(file))
    stop('you did not run rms::stanCompile to compile Stan code')
  mod <- readRDS(file)

  opt <- NULL
  if(method != 'sampling') {
    otime <- system.time(
      g <- rstan::optimizing(mod, data=d, ...))
    if(g$return_code != 0)
      warning(paste('optimizing did not work; return code', g$return_code))
    parm <- g$par
    nam <- names(parm)
    th  <- nam[grep('theta\\[', nam)]
    al  <- nam[grep('alpha\\[', nam)]
    be  <- nam[grep('beta\\[', nam)]
    alphas <- parm[al]
    betas  <- parm[be]
    thetas <- parm[th]
    names(alphas) <- if(nrp == 1) 'Intercept' else paste0('y>=', ylev[-1])
    alphas <- alphas - sum(betas * xbar)
    names(betas)  <- names(thetas) <- atr$colnames
    opt <- list(coefficients=c(alphas, betas), theta=thetas,
              sigmag=parm['sigmag'], deviance=-2 * g$value,
              return_code=g$return_code, hessian=g$hessian,
              executionTime=otime)
    if(method == 'optimizing') return(opt)
    }

	g <- rstan::sampling(mod, pars='sigmaw', include=FALSE,
                       data=d, iter=iter, chains=chains, refresh=refresh, ...)
  ## g <- rstan::vb(mod, data=d)   # bombed
	nam <- names(g)
	al  <- nam[grep('alpha\\[', nam)]
	be  <- nam[grep('beta\\[', nam)]
  ga  <- nam[grep('gamma\\[', nam)]

  draws  <- as.matrix(g)
	alphas <- draws[, al, drop=FALSE]
	betas  <- draws[, be, drop=FALSE]

  sigmags <- NULL; clparm <- character(0); gammas <- NULL
  if(length(cluster)) {
    sigmags <- draws[, 'sigmag']
    clparm  <- 'sigmag'
    cle     <- draws[, ga, drop=FALSE]
    gammas  <- apply(cle, 2, median)
  }

  epsmed <- rhos <- NULL
  if(length(time)) {
    rhos    <- draws[, 'rho']
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
  
	diagnostics <-
		rstan::summary(g, pars=c(al, be, clparm),
                   probs=NULL)$summary[,c('n_eff', 'Rhat')]
	# Back-scale to original data scale
	alphacorr <- rowSums(sweep(betas, 2, xbar, '*')) # was xbar/xsd
	alphas    <- sweep(alphas, 1, alphacorr, '-')
	# betas     <- sweep(betas, 2, xsd, '/')
	colnames(alphas) <- if(nrp == 1) 'Intercept' else paste0('y>=', ylev[-1])
	colnames(betas)  <- atr$colnames
	draws            <- cbind(alphas, betas)

  param <- rbind(mean=colMeans(draws), median=apply(draws, 2, median))
  if(method != 'sampling') {
    param <- rbind(mode=opt$coefficients, param)
    opt$coefficients <- NULL
  }

  Loo <- if(loo) rstan::loo(g)

  freq <- table(Y, dnn=yname)
    
	res <- list(call=call, 
							draws=draws, sigmags=sigmags, rhos=rhos,
              gammas=gammas, eps=epsmed,
              param=param, N=n, p=p, yname=yname, ylevels=ylev, freq=freq,
						  alphas=al, betas=be,
						  xbar=xbar, Design=atr, scale.pred=c('log odds', 'Odds Ratio'),
							terms=Terms, assign=ass, na.action=atrx$na.action, fail=FALSE,
							non.slopes=nrp, interceptRef=kmid, sformula=sformula,
							x=if(x) X, y=if(y) Y, loo=Loo,
              clusterInfo=if(length(cluster))
                list(cluster=if(x) cluster else NULL, n=Nc, name=clustername),
              timeInfo=if(length(time))
                         list(time=if(x) time else NULL, n=Nt, nobs=Ntobs,
                              name=timename),
							rstan=g, opt=opt, diagnostics=diagnostics,
              iter=iter, chains=chains)
	class(res) <- c('blrm', 'rmsb', 'rms')
	res
}

##' Compute Indexes of Predictive Accuracy and Their Uncertainties
##'
##' For a binary or ordinal logistic regression fit from \code{blrm}, computes several indexes of predictive accuracy along with credible intervals for them.  Optionally plots their posterior densities.
##' @title blrmStats
##' @param fit an object produced by \code{blrm}
##' @param ns number of posterior draws to use in the calculations (default is 400)
##' @param cint credible interval probability (default is 0.95)
##' @param pl set to \code{TRUE} to plot the posterior densities using base graphics
##' @return list of class \code{'blrmStats'} whose most important element is \code{Stats}.  The indexes computed are defined below, with gp, B, EV, and vp computed using the intercept corresponding to the median value of Y.  See \url{http://fharrell.com/post/add-value} for more information.
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
blrmStats <- function(fit, ns=400, cint=0.95, pl=FALSE) {
  X <- fit$x
  y <- fit$y
  if(length(X) == 0 | length(y) == 0)
    stop('must have specified x=TRUE, y=TRUE to blrm')
  y <- as.integer(y) - 1
  nrp    <- fit$non.slopes
  kmid   <- fit$interceptRef  # intercept close to median
  alphas <- fit$draws[, kmid]
  ## Save betas as a matrix
  s      <- fit$draws[, - (1 : nrp), drop=FALSE]  # omit intercepts
  ndraws <- nrow(s)
  ns     <- min(ndraws, ns)
  if(ns < ndraws) {
    j <- sample(1 : ndraws, ns, replace=FALSE)
    s <- s[j,, drop=FALSE]
    alphas <- alphas[j]
    }
  ylev  <- fit$ylevels
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


  for(i in 1 : ns) {
    beta  <- s[i, ]
    alpha <- alphas[i]
    lp   <- alpha + X %*% beta
    prob <- plogis(lp)
    d    <- dxy(lp, y)
    C    <- (d + 1.) / 2.
    st <- c(d, C, GiniMd(lp), GiniMd(prob),
            brier(prob, y > kmid - 1), 
            br2(prob), var(lp), var(prob))
    stats[i,] <- st
  }

  sbar <- colMeans(stats)
  se   <- apply(stats, 2, sd)
  alpha <- 1 - cint
  cred <- apply(stats, 2, quantile,
                probs=c(alpha / 2, 1 - alpha / 2, .05, .95))
  sym  <- (cred[4, ] - sbar) / (sbar - cred[3, ])
  
  text <- paste0(round(sbar, 3), ' [', round(cred[1,], 3), ', ',
                 round(cred[2,], 3), ']')
  Stats <- rbind(Mean=sbar, SE=se, ci=cred[1:2,], Symmetry=sym)
  statnames <- colnames(Stats)
  names(text) <- statnames

  if(pl) {
    par(mfrow=c(4,2), mar=c(3, 2, 0.5, 0.5), mgp=c(1.75, .55, 0))
    for(w in setdiff(statnames, 'C')) {
      den <- density(stats[, w])
      plot(den, xlab=w, ylab='', type='l', main='')
      ref <- c(sbar[w], cred[1, w], cred[2, w])
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
##' Prints main results from \code{blrm} along with indexes and predictive accuracy and their credible intervals computed from \code{blrmStats}.
##' @title print.blrm
##' @param x object created by \code{blrm}
##' @param dec number of digits to print to the right of the decimal
##' @param coefs specify \code{FALSE} to suppress printing parameter estimates, and in integer k to print only the first k
##' @param cint credible interval probability for summary indexes
##' @param ns number of random samples of the posterior draws for use in computing credible intervals for accuracy indexes
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
print.blrm <- function(x, dec=4, coefs=TRUE, cint=0.95, ns=400,
                      title='Bayesian Logistic Regression Model', ...) {
  latex <- prType() == 'latex'
  
  z <- list()
  k <- 0
  
  if(length(x$freq) > 3) {
    k <- k + 1
    z[[k]] <- list(type='print', list(x$freq),
                   title='Frequencies of Responses')
  }
  if(length(x$na.action)) {
    k <- k + 1
    z[[k]] <- list(type=paste('naprint',class(x$na.action),sep='.'),
                   list(x$na.action))
  }

  alp <- (1. - cint) / 2.
  qro <- function(x) {
    r <- round(quantile(x, c(alp, 0.5, 1. - alp)), 4)
    paste0(r[2], ' [', r[1], ', ', r[3], ']')
    }
  ci  <- x$clusterInfo
  sigmasum <- NULL
  if(length(ci)) sigmasum <- qro(x$sigmags)

  ti  <- x$timeInfo
  rhosum <- if(length(ti)) qro(x$rhos)
 
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
                      p               = x$p,
                      'Cluster on'    = ci$name,
                      Clusters        = ci$n,
                      'sigma gamma'   = sigmasum,
                      'Time variable' = ti$name,
                      'Max time'      = ti$n,
                      'Obs times'     = ti$nobs,
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
                   list(bayes=print.rmsb(x, cint=cint, pr=FALSE)))
  }
  
  prModFit(x, title=title, z, digits=dec, coefs=coefs, ...)
}

## ??
## Code for mean below needs to be finished adapting to Bayes
  
##' Make predictions from a \code{blrm} fit
##'
##' Predict method for \code{blrm} objects
##' @title predict.blrm
##' @param object,...,type,se.fit,codes see [predict.lrm] 
##' @param posterior.summary set to \code{'median'} or \code{'mode'} to use posterior median/mode instead of mean
##' @return a data frame,  matrix, or vector
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   predict(f, newdata, posterior.summary='median')
##' }
##' @seealso [predict.lrm]
##' @author Frank Harrell
##' @md
predict.blrm <- function(object, ..., 
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		  "terms", "cterms", "ccterms", "adjto", "adjto.data.frame",
      "model.frame"),
		se.fit=FALSE, codes=FALSE,
    posterior.summary=c('mean', 'median', 'mode')) {
  
  type           <- match.arg(type)
  posterior.summary <- match.arg(posterior.summary)

  if(type %in% c('x', 'data.frame', 'terms', 'cterms', 'ccterms',
                 'adjto', 'adjto.data.frame', 'model.frame'))
    return(predictrms(object,...,type=type,
                      posterior.summary=posterior.summary))

  if(type != 'x') stop('types other than x not yet implemented')
  X <- predictrms(object, ..., type='x', posterior.summary=posterior.summary)
  lp <- X %*% object$draws
  
  xb <- predictrms(object, ..., type="lp", se.fit=FALSE,
                   posterior.summary=posterior.summary)
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
