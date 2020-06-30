##' Print Stan Diagnostics
##'
##' Retrieves the effect samples sizes and Rhats computed after a fitting function ran \code{rstan}, and prepares it for printing.  If the fit was created by \code{stackImpute}, the diagnostics for all imputations are printed (separately).
##' @title stanDx
##' @param object an object created by an \code{rms} package Bayesian fitting function such as \code{blrm} or \code{stackMI}
##' @return matrix suitable for printing
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   stanDx(f)
##' }
##' @author Frank Harrell
stanDx <- function(object) {
  draws    <- object$draws
  n.impute <- object$n.impute
  if(! length(n.impute)) n.impute <- 1
  d <- object$diagnostics
  if(n.impute == 1 && (length(names(d)) == 1 && names(d) == 'pars')) {
    cat('Diagnostics failed trying to summarize', d$pars, '\n')
    return(invisible())
  }
  if(n.impute > 1) cat('Diagnostics for each of', n.impute, 'imputations\n\n')
	cat('Iterations:', object$iter, 'on each of', object$chains, 'chains, with',
			nrow(draws) / n.impute, 'posterior distribution samples saved\n\n')
	cat('For each parameter, n_eff is a crude measure of effective sample size',
			'and Rhat is the potential scale reduction factor on split chains',
			'(at convergence, Rhat=1)\n', sep='\n')
  if(n.impute == 1) {
    d[, 'n_eff']  <- round(d[, 'n_eff'])
    d[, 'Rhat']   <- round(d[, 'Rhat'], 3)
    rownames(d)   <- colnames(cbind(draws, object$omega))
    return(d)
  } else {
    D <- NULL
    for(i in 1 : n.impute) {
      dx <- d[[i]]
      dx[, 'n_eff']  <- round(dx[, 'n_eff'])
      dx[, 'Rhat' ]  <- round(dx[, 'Rhat'], 3)
      nams <- colnames(cbind(draws, object$omega))
      rownames(dx)   <- paste0(paste0('Imputation ', i, ': '), nams)
      D <- rbind(D, dx)
      }
  }
  D
}

##' Get Stan Output
##'
##' Extracts the object created by \code{rstan::sampling} so that standard Stan diagnostics can be run from it
##' @title stanGet
##' @param object an objected created by an \code{rms} package Bayesian fitting function
##' @return the object created by \code{rstan::sampling}
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   s <- stanGet(f)
##' }
##' @author Frank Harrell
stanGet <- function(object) object$rstan

##' Extract Bayesian Summary of Coefficients
##'
##' Computes either the posterior mean (default), posterior median, or posterior mode of the parameters in an \code{rms} Bayesian regression model
##' @title coef.rmsb
##' @param object an object created by an \code{rms} package Bayesian fitting function
##' @param stat name of measure of posterior distribution central tendency to compute
##' @param ... ignored
##' @return a vector of intercepts and regresion coefficients
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   coef(f, stat='mode')
##' }
##' @author Frank Harrell
coef.rmsb <- function(object, stat=c('mean', 'median', 'mode'), ...) {
	stat <- match.arg(stat)
	# pmode <- function(x) {
	# 	dens <- density(x)
	# 	dens$x[which.max(dens$y)[1]]
	# }
  getParamCoef(object, stat)
}

##' Variance-Covariance Matrix
##'
##' Computes the variance-covariance matrix from the posterior draws by compute the sample covariance matrix of the draws
##' @title vcov.rmsb
##' @param object an object produced by an \code{rms} package Bayesian fitting function
##' @param regcoef.only set to \code{FALSE} to also include non-regression coefficients such as shape/scale parameters
##' @param intercepts set to \code{'all'} to include all intercepts (the default), \code{'none'} to exclude them all, or a vector of integers to get selected intercepts
##' @param ... ignored
##' @return matrix
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   v <- vcov(f)
##' }
##' @seealso [vcov.rms]
##' @author Frank Harrell
##' @md
vcov.rmsb <- function(object, regcoef.only=TRUE,
                      intercepts='all', ...) {
  
  ## Later will have to handle non-coefficient parameters

  if(length(intercepts) == 1 && is.character(intercepts) &&
     intercepts %nin% c('all', 'none'))
    stop('if character, intercepts must be "all" or "none"')

  draws <- object$draws
  
  if(! length(intercepts) ||
     (length(intercepts) == 1) && intercepts == 'all')
    return(var(draws))

  ns <- num.intercepts(object)
  p  <- ncol(draws)
  nx <- p - ns
  if(intercepts == 'none') intercepts <- integer(0)
  i <- if(nx == 0) intercepts else c(intercepts, (ns + 1) : p)
  var(draws[, i, drop=FALSE])
}

##' Basic Print for Bayesian Parameter Summary
##'
##' For a Bayesian regression fit prints the posterior mean, median, SE, highest posterior density interval, and symmetry coefficient from the posterior draws.  For a given parameter, the symmetry measure is computed using the \code{distSym} function.
##' @title print.rmsb
##' @param x an object created by an \code{rms} Bayesian fitting function
##' @param prob HPD interval coverage probability (default is 0.95)
##' @param dec amount of rounding (digits to the right of the decimal)
##' @param intercepts set to \code{FALSE} to not print intercepts
##' @param pr set to \code{FALSE} to return an unrounded matrix and not print
##' @param ... ignored
##' @return matrix (rounded if \code{pr=TRUE})
##' @examples
##' \dontrun{
##'   f <- blrm(...)
##'   print.rmsb(f)
##' }
##' @author Frank Harrell
print.rmsb <- function(x, prob=0.95, dec=4, intercepts=TRUE, pr=TRUE, ...) {
  nrp   <- num.intercepts(x)
	s     <- x$draws
  param <- t(x$param)
  if(! intercepts && nrp > 0) {
    s     <- s[,   -(1 : nrp),  drop=FALSE]
    param <- param[-(1 : nrp),, drop=FALSE]
    }
  means <- param[, 'mean']
  colnames(param) <- upFirst(colnames(param))

	se  <- sqrt(diag(var(s)))
	hpd <- apply(s, 2, HPDint, prob=prob)
  P   <- apply(s, 2, function(u) mean(u > 0))
	sym <- apply(s, 2, distSym)
	w <- cbind(param, SE=se, Lower=hpd[1,], Upper=hpd[2,], P, Symmetry=sym)
  rownames(w) <- names(means)
  if(! pr) return(w)
	cat(nrow(s), 'draws from the posterior distribution\n\n')
	round(w, dec)
}

##' Compile Stan Code
##'
##' Retrieves Stan code files from the github repository, compiles then with \code{rstan::stan_model}, and stores the compiled code in a central placed defined by the user with e.g. \code{options(stancompiled='~/R/stan')}.
##' stanCompile
##' @param mods character vector of base names of \code{.stan} model files to compile.  Default is to compile all those currently in use in \code{rms}
##' @param repo URL to online source file base
##' @examples
##'   \dontrun{
##'   options(stancompiled='~/R/stan')    # need this always
##'   stanCompile()    # do this once per computer to compile centrally
##'   }
##' @author Frank Harrell
stanCompile <-
  function(mods=c('lrm', 'lrmc', 'lrmppo', 'lrmcppo', 'lrmcar1'),
           repo='https://raw.githubusercontent.com/harrelfe/stan/master') {
  requireNamespace('rstan', quietly=TRUE)
  options(auto_write = FALSE)
  stanloc <- .Options$stancompiled
  if(! length(stanloc)) stop('options(stancompiled=) not defined')

  cat('Compiling', length(mods), 'programs to', stanloc, '\n')
  for(m in mods) {
    cat('Compiling', m, '\n')
    f <- paste0(repo, '/', m, '.stan')
    w <- readLines(f)
    k <- rstan::stan_model(model_code = w)
    mod <- paste0(stanloc, '/', m, '.rds')
    saveRDS(k, mod)
  }
  invisible()
  }

##' Plot Posterior Densities and Summaries
##'
##' For an \code{rms} Bayesian fit object, plots posterior densities for selected parameters along with posterior mode, mean, median, and highest posterior density interval.  If the fit was produced by \code{stackMI} the density represents the distribution after stacking the posterior draws over imputations, and the per-imputation density is also drawn as pale curves.  If exactly two parameters are being plotted and \code{bivar=TRUE}, hightest bivariate posterior density contrours are plotted instead, for a variety of \code{prob} values including the one specified, using 
##' @title plot.rmsb
##' @param x an \code{rms} Bayesian fit object
##' @param which names of parameters to plot, defaulting to all non-intercepts. Can instead be a vector of integers.
##' @param nrow number of rows of plots
##' @param ncol number of columns of plots
##' @param prob probability for HPD interval
##' @param bivar set to \code{TRUE} to plot bivariate density contours instead of univariate results (ignored if the number of parameters plotted is not exactly two)
##' @param bivarmethod passed as \code{method} argument to \code{pdensityContour}
##' @param ... passed to \code{pdensityContour}
##' @return \code{ggplot2} object
##' @author Frank Harrell
plot.rmsb <- function(x, which=NULL, nrow=NULL, ncol=NULL, prob=0.95,
                      bivar=FALSE, bivarmethod=c('ellipse', 'kernel'), ...) {

  bivarmethod <- match.arg(bivarmethod)
  
  nrp      <- num.intercepts(x)
  draws    <- x$draws

  n.impute <- x$n.impute
  if(! length(n.impute)) n.impute <- 1

  omega <- x$omega
  est   <- x$param
  if(length(omega)) {
    draws <- cbind(draws, omega)
    g <- function(x) c(mode=NA, mean=mean(x), median=median(x))
    osummary <- apply(omega, 2, g)
    osummary <- osummary[rownames(est),, drop=FALSE]
    est      <- cbind(est, osummary)
    }
  nd    <- nrow(draws)
  nam   <- colnames(draws)
  if(! length(which)) which <- if(nrp == 0) nam else nam[-(1 : nrp)]
  if(! is.character(which)) which <- nam[which]
  draws <- draws[, which, drop=FALSE]

  if(length(which) == 2 && bivar) {
    g    <- pdensityContour(draws[, 1], draws[, 2], prob=prob, pl=TRUE,
                            method=bivarmethod, ...)
    cn   <- colnames(draws)
    g    <- g + xlab(cn[1]) +  ylab(cn[2])
    return(g)
    }
  
  hpd   <- apply(draws, 2, HPDint, prob=prob)
  
  draws  <- as.vector(draws)
  param  <- factor(rep(which, each=nd), which)
  imputation <- rep(rep(1 : n.impute, each = nd / n.impute), length(which))
  imputation <- paste('Imputation', imputation)

  est    <- est[, which, drop=FALSE]
  est    <- rbind(est, hpd)
  ne     <- nrow(est)
  stat   <- rownames(est)
  stat   <- ifelse(stat %in% c('Lower', 'Upper'),
                   paste(prob, 'HPDI'), stat)
  est    <- as.vector(est)
  eparam <- factor(rep(which, each=ne), which)
  stat   <- rep(stat, length(which))

  sz <- if(n.impute > 1) 0.2 else 1.4
  d  <- data.frame(param, draws,
                   imputation = if(n.impute > 1) imputation else 'Density',
                   sz         = sz)
  if(n.impute > 1) {
    d2 <- d
    d2$imputation <- 'Stacked'
    d2$sz <- 1.4
    d <- rbind(d, d2)
    implev <- c('Stacked', paste('Imputation', 1:n.impute))
    d$imputation <- factor(d$imputation, implev, implev)  # ignored ??
  }

  de <- data.frame(param=eparam, est, stat)
  g <- ggplot(d, aes(x=draws, color=imputation, size=I(sz))) + geom_density() +
         geom_vline(data=de, aes(xintercept=est, color=stat, alpha=I(0.4))) +
         facet_wrap(~ param, scales='free', nrow=nrow, ncol=ncol) +
         guides(color=guide_legend(title='')) +
         xlab('') + ylab('')
  g
}

##' Diagnostic Trace Plots
##'
##' For an \code{rms} Bayesian fit object, uses by default the stored posterior draws to check convergence properties of posterior sampling.  If instead \code{rstan=TRUE}, calls the \code{rstan} \code{traceplot} function on the \code{rstan} object inside the \code{rmsb} object, to check properties of posterior sampling.  If \code{rstan=TRUE} and the \code{rstan} object has been removed and \code{previous=TRUE}, attempts to find an already existing plot created by a previous run of the \code{knitr} chunk, assuming it was the \code{plotno} numbered plot of the chunk.
##' @title stanDxplot
##' @param x an \code{rms} Bayesian fit object
##' @param which names of parameters to plot, defaulting to all non-intercepts.  When \code{rstan=FALSE} these are the friendly \code{rms} names, otherwise they are the \code{rstan} parameter names.  If the model fit was run through \code{stackMI} for multiple imputation, the number of traces is multiplied by the number of imputations.  Set to \code{'ALL'} to plot all parameters.
##' @param rstan set to \code{TRUE} to use \code{rstan::traceplot} on a (presumed) stored \code{rstan} object in \code{x}, otherwise only real iterations are plotted and parameter values are shown as points instead of lines, with chains separated
##' @param previous see details
##' @param plotno see details
##' @param rev set to \code{TRUE} to reverse direction for faceting chains
##' @param stripsize specifies size of chain facet label text, default is 8
##' @param ... passed to \code{rstan::traceplot}
##' @return \code{ggplot2} object if \code{rstan} object was in \code{x}
##' @author Frank Harrell
stanDxplot <- function(x, which=NULL, rstan=FALSE, previous=TRUE,
                       plotno=1, rev=FALSE, stripsize=8, ...) {
  if(! rstan) {
    draws <- x$draws
    if(length(which) == 1 && which == 'ALL') which <- colnames(draws)
    if(! length(which)) {
      nrp <- x$non.slopes
      which <- if(nrp > 0) colnames(draws)[-(1 : nrp)] else colnames(draws)
    }
    draws   <- cbind(draws, x$omega)
    draws   <- draws[, which, drop=FALSE]
    n.impute <- x$n.impute
    if(! length(n.impute)) n.impute <- 1
    nchains <- x$chains * n.impute
    ndraws  <- nrow(draws)
    chain   <- rep(rep(1 : nchains, each = ndraws / nchains), length(which))
    chain   <- paste('Chain', chain)
    draws   <- as.vector(draws)
    param   <- rep(which, each=ndraws)
    iter    <- rep(rep(1 : (ndraws / nchains), nchains), length(which))
    d       <- data.frame(chain, iter, param, draws)
    g <- ggplot(d, aes(x=iter, y=draws)) +
      geom_point(size=I(0.03), alpha=I(0.3)) +
      xlab('Post Burn-in Iteration') + ylab('Parameter Value') +
      theme(strip.text = element_text(size = stripsize))

    if(rev)  g <- g + facet_grid(chain ~ param, scales='free')
    else     g <- g + facet_grid(param ~ chain, scales='free_y')
    return(g)
    }

  if(! length(which)) which <- x$betas
  s <- x$rstan
  if(length(s)) return(rstan::traceplot(s, pars=which, ...))

  if(! previous)
    stop('fit did not have rstan information and you specified previous=FALSE')
  ## Assume that the plot was generated in a previous run before
  get   <- knitr::opts_current$get
  path  <- get('fig.path')
  cname <- get('label')
  dev   <- get('dev')
  if(! length(cname))
    stop('with no rstan component in fit and previous=TRUE you must call stanDxplot inside a knitr Rmarkdown chunk')
  file <- paste0(path, cname, '-', plotno, '.', dev)
  if(! file.exists(file))
    stop(paste('with no rstan component in fit, file', file, 'must exist when previous=TRUE'))
  knitr::include_graphics(file)
  }


##' Function Generator for Posterior Probabilities of Assertions
##'
##' From a Bayesian fit object such as that from \code{blrm} generates an R function for evaluating the probability that an assertion is true.  The probability, within simulation error, is the proportion of times the assertion is true over the posterior draws.  If the assertion does not evaluate to a logical or 0/1 quantity, it is taken as a continuous derived parameter and the vector of draws for that parameter is returned and can be passed to the \code{PostF} plot method.  \code{PostF} can also be used on objects created by \code{contrast.rms}
##' @title PostF
##' @param fit a Bayesian fit or \code{contrast.rms} object
##' @param name specifies whether assertions will refer to shortened parameter names (the default) or original names.  Shorted names are of the form \code{a1, ..., ak} where \code{k} is the number of intercepts in the model, and \code{b1, ..., bp} where \code{p} is the number of non-intercepts.  When using original names that are not legal R variable names, you must enclose them in backticks.  For \code{contrast} objects, \code{name} is ignored and you must use contrast names.  The \code{cnames} argument to \code{contrast.rms} is handy for assigning your own names.
##' @param pr set to \code{TRUE} to have a table of short names and original names printed when \code{name='short'}.  For \code{contrasts} the contrast names are printed if \code{pr=TRUE}.
##' @return an R function
##' @examples
##' \dontrun{
##'   f <- blrm(y ~ age + sex)
##'   P <- PostF(f)
##'   P(b2 > 0)     # Model is a1 + b1*age + b2*(sex == 'male')
##'   P(b1 < 0 & b2 > 0)   # Post prob of a compound assertion
##'   # To compute probabilities using original parameter names:
##'   P <- PostF(f, name='orig')
##'   P(age < 0)    # Post prob of negative age effect
##'   P(`sex=male` > 0)
##'   f <- blrm(y ~ sex + pol(age, 2))
##'   P <- PostF(f)
##'   # Compute and plot posterior density of the vertex of the
##'   # quadratic age effect
##'   plot(P(-b2 / (2 * b3)))
##'
##'   # The following would be useful in age and sex interacted
##'   k <- contrast(f, list(age=c(30, 50), sex='male'),
##'                    list(age=c(30, 50), sex='female'),
##'                 cnames=c('age 30 M-F', 'age 50 M-F'))
##'   P <- PostF(k)
##'   P(`age 30 M-F` > 0 & `age 50 M-F` > 0)
##' ##' }
##' @author Frank Harrell
PostF <- function(fit, name=c('short', 'orig'), pr=FALSE) {
  name       <- match.arg(name)
  alphas     <- fit$alphas
  betas      <- fit$betas
  draws      <- fit$draws
  # See if fit is really an object created by contrast.rms
  iscon      <- length(fit$cdraws) > 0
  if(iscon) draws <- fit$cdraws
  if(iscon && pr) cat('Contrast names:', paste(colnames(draws), collapse=', '),
                      '\n')
  
  orig.names <- colnames(draws)
  if(! iscon && (name == 'short')) {
    nrp <- num.intercepts(fit)
    rp  <- length(orig.names) - nrp
    nam <- c(if(nrp > 0) paste0('a', 1 : nrp),
             paste0('b', 1 : rp))
    if(pr) {
      w <- cbind('Original Name' = orig.names,
                 'Short Name'    = nam)
      rownames(w) <- rep('', nrow(w))
      print(w, quote=FALSE)
      }
    colnames(draws) <- nam
    }
  f <- function(assert, label, draws) {
    if(! length(label)) label <- as.character(sys.call()[2])
    w <- eval(substitute(assert), draws)
    if(length(unique(w)) < 3) return(mean(w))
    structure(w, class='PostF', label=label)
    }

  # Convert draws to data frame so eval() will work
  formals(f) <- list(assert=NULL, label=NULL, draws=as.data.frame(draws))
  f
}

##' Plot Posterior Density of \code{PostF}
##'
##' Computes highest posterior density and posterior mean and median as vertical lines, and plots these on the density function.  You can transform the posterior draws while plotting.
##' @title plot.PostF
##' @param x result of running a function created by \code{PostF}
##' @param ... other results created by such functions
##' @param cint interval probability
##' @param label x-axis label if not the expression originally evaluated.  When more than one result is plotted, \code{label} is a vector of character strings, one for each result.
##' @param type when plotting more than one result specifies whether to make one plot distinguishing results by line type, or whether to make separate panels
##' @param ltitle used of \code{type='linetype'} to specify name of legend for the line types
##' @return \code{ggplot2} object
##' @author Frank Harrell
plot.PostF <- function(x, ..., cint=0.95, label=NULL,
                       type=c('linetype', 'facet'), ltitle='') {
  clab  <- paste(cint, 'HPD\nInterval')
  type  <- match.arg(type)
  d <- list(...)
  if(! length(d)) {
    if(! length(label)) label <- attr(x, 'label')
    x     <- as.vector(x)
    hpd   <- HPDint(x, cint)
    de    <- data.frame(est=c(mean(x), median(x), hpd),
                        stat=c('Mean', 'Median', rep(clab, 2)))
    g <- ggplot(mapping=aes(x=x)) + geom_density() +
      geom_vline(data=de, aes(xintercept=est, color=stat, alpha=I(0.4))) +
      guides(color=guide_legend(title='')) +
      xlab(label) + ylab('')
    return(g)
  }

  d <- c(list(x), d)
  if(! length(label)) label <- sapply(d, function(x) attr(x, 'label'))
  if(length(label) != length(d)) stop('label has incorrect length')
  parm <- rep(label, sapply(d, length))
  X    <- unlist(d)
  de   <- NULL
  for(i in 1 : length(d)) {
    x   <- as.vector(d[[i]])
    hpd <- HPDint(x, cint)
    w   <- data.frame(est=c(mean(x), median(x), hpd),
                      stat=c('Mean', 'Median', rep(clab, 2)),
                      parm=label[i])
    de  <- rbind(de, w)
  }
  switch(type,
         linetype = {
           ggplot(mapping=aes(x=X, linetype=parm)) + geom_density() +
             geom_vline(data=de, aes(xintercept=est, color=stat, alpha=I(0.4),
                                     linetype=parm)) +
             guides(color=guide_legend(title=''),
                    linetype=guide_legend(title=ltitle)) +
             xlab('') + ylab('')
         },
         facet = {
           w <- data.frame(X, parm)
           ggplot(w, mapping=aes(x=X)) + geom_density() + facet_wrap(~ parm) +
             geom_vline(data=de, aes(xintercept=est, color=stat, alpha=I(0.4))) +
             guides(color=guide_legend(title='')) +
             xlab('') + ylab('')
           })
}


## Test this approach:
#f <- function(assert, draws) eval(substitute(assert), draws)
#formals(f) <- list(assert=NULL, draws=list(a=1:5, b1=2:6, b3=3:7))
#f(a); f(b1); f(b2)   # with(draws, assert) did not work

# Help info is in rmsMisc.Rd since this is a non-exported function
getParamCoef <- function(fit, posterior.summary=c('mean', 'median', 'mode')) {
  posterior.summary <- match.arg(posterior.summary)
  param <- fit$param
  if(posterior.summary == 'mode' && 'mode' %nin% rownames(param))
    stop('posterior mode not included in model fit')
  param[posterior.summary, ]
}


##' Compare Bayesian Model Fits
##'
##' Uses \code{loo::loo_model_weights} to compare a series of models such as those created with \code{blrm}
##' @title compareBmods
##' @param ... a series of model fits
##' @param method see [loo::loo_model_weights]
##' @param r_eff_list see [loo::loo_model_weights]
##' @return a \code{loo::loo_model_weights} object
##' @author Frank Harrell
##' @md
compareBmods <- function(..., method='stacking', r_eff_list=NULL) {
  fits <- list(...)
  lo   <- lapply(fits, function(x) x$loo)
  loo::loo_model_weights(lo, method=method, r_eff_list=r_eff_list)
  }

##' Highest Posterior Density Interval
##'
##' Adapts code from \code{coda::HPDinterval} to compute a highest posterior density interval from posterior samples for a single parameter.  Quoting from the \code{coda} help file, for each parameter the interval is constructed from the empirical cdf of the sample as the shortest interval  for  which  the  difference  in  the  ecdf  values  of  the  endpoints  is  the  nominal  probability.  Assuming that the distribution is not severely multimodal, this is the HPD interval.
##' @title HPDint 
##' @param x a vector of posterior draws
##' @param prob desired probability coverage
##' @return a 2-vector with elements \code{Lower} and \code{Upper}
##' @author Douglas Bates and Frank Harrell
HPDint <- function(x, prob = 0.95) {
  x     <- sort(x)
  nsamp <- length(x)
  gap   <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init  <- 1:(nsamp - gap)
  inds  <- which.min(x[init + gap] - x[init])
  c(Lower =  x[inds], Upper = x[inds + gap])
}

##' Distribution Symmetry Measure
##'
##' From a sample from a distribution computes a symmetry measure.  By default it is the gap between the mean and the 0.95 quantile divided by the gap between the 0.05 quantile and the mean.
##' @title distSym
##' @param x a numeric vector representing a sample from a continuous distribution
##' @param prob quantile interval coverage
##' @param na.rm set to \code{TRUE} to remove \code{NA}s before proceeding.
##' @return a scalar with a value of 1.0 indicating symmetry
##' @author Frank Harrell
distSym <- function(x, prob=0.9, na.rm=FALSE) {
  if(na.rm) x <- x[! is.na(x)]
  a <- (1. - prob) / 2.
  w <- quantile(x, probs=c(a / 2., 1. - a / 2.))
  xbar <- mean(x)
  (w[2] - xbar) / (xbar - w[1])
}

##' Bivariate Posterior Contour
##'
##' Computes coordinates of a highest density contour containing a given probability volume given a sample from a continuous bivariate distribution, and optionally plots.  The default method assumes an elliptical shape, but one can optionally use a kernel density estimator.
##' Code adapted from \code{embbook::HPDregionplot}.  See \url{http://www.sumsar.net/blog/2014/11/how-to-summarize-a-2d-posterior-using-a-highest-density-ellipse}.
##' @title pdensity2d
##' @param x a numeric vector 
##' @param y a numeric vector the same length of x
##' @param prob main probability coverage (the only one for \code{method='ellipse'})
##' @param otherprob vector of other probability coverages for \code{method='kernel'}
##' @param method defaults to \code{'ellipse'}, can be set to \code{'kernel'}
##' @param h vector of bandwidths for x and y.  See \code{MASS::kde2d}.
##' @param n number of grid points in each direction, defaulting to normal reference bandwidth (see \code{bandwidth.nrd}).
##' @param pl set to \code{TRUE} to plot contours
##' @return a 2-column matrix with x and y coordinates unless \code{pl=TRUE} in which case a \code{ggplot2} graphic is returned
##' @author Ben Bolker and Frank Harrell
pdensityContour <-
  function(x, y, method=c('ellipse', 'kernel'),
           prob=0.95, otherprob=c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9),
           h=c(1.3 * MASS::bandwidth.nrd(x),
               1.3 * MASS::bandwidth.nrd(y)),
           n=70, pl=FALSE) {

method <- match.arg(method)

rho    <- cor(x, y, method='spearman')
rholab <- bquote(paste('Spearman ', rho == .(round(rho, 2))))
    
if(rho > 0.999) {
  r    <- range(x)
  xout <- seq(r[1], r[2], length=150)
  lo   <- lowess(x, y)
  d    <- as.data.frame(approx(lo, xout=xout))
  if(pl) {
    g <- ggplot(d, aes(x=x, y=y)) + geom_line() + labs(caption=rholab)
    return(g)
  }
  return(d)
}
    
if(method == 'ellipse') {
  xy <- cbind(x, y)
  f <- MASS::cov.mve(xy, quantile.used=round(nrow(xy) * prob))
  points_in_ellipse <- xy[f$best, ]
  boundary <- predict(cluster::ellipsoidhull(points_in_ellipse))
  d <- data.frame(x=boundary[, 1], y=boundary[, 2])
  if(pl) {
 
    g <- ggplot(d, aes(x=x, y=y)) + geom_path() + labs(caption=rholab)
    return(g)
  }
  return(d)
}

  prob <- unique(sort(c(prob, otherprob)))
    
  f <- MASS::kde2d(x, y, n=n, h=h)
  x <- f$x
  y <- f$y
  z <- f$z
  
  dx <- diff(x[1:2])
  dy <- diff(y[1:2])
  sz <- sort(z)
  c1 <- cumsum(sz) * dx * dy
  ## trying to find level containing prob of volume ...
  levels <- sapply(prob, function(p) approx(c1, sz, xout = 1. - p, rule=2)$y)
  X <- Y <- Prob <- NULL
  for(i in 1 : length(prob)) {
    w <- contourLines(x, y, z, levels=levels[i])
    X <- c(X, w[[1]]$x)
    Y <- c(Y, w[[1]]$y)
    Prob <- c(Prob, rep(prob[i], length(w[[1]]$x)))
  }
    d <- data.frame(x=X, y=Y, Probability = factor(Prob))
  if(pl)
    ggplot(d, aes(x=x, y=y, col=Probability)) + geom_path() +
      scale_color_brewer(direction = -1) + labs(caption=rholab)
  else d
  }
utils::globalVariables('Probability')     # why in the world needed?

##' QR Decomposition Preserving Selected Columns
##'
##' Runs a matrix through the QR decomposition and returns the transformed matrix and the forward and inverse transforming matrices \code{R, Rinv}.  If columns of the input matrix \code{X} are centered the QR transformed matrix will be orthogonal.  This is helpful in understanding the transformation and in scaling prior distributions on the transformed scale.  \code{not} can be specified to keep selected columns as-is.  \code{cornerQr} leaves the last column of \code{X} alone (possibly after centering).  When \code{not} is specified, the square transforming matrices have appropriate identity submatrices inserted so that recreation of original \code{X} is automatic.
##' @title selectedQr 
##' @param X a numeric matrix
##' @param not an integer vector specifying which columns of \code{X} are to be kept with their original values
##' @param corner set to \code{FALSE} to not treat the last column specially.  You may not specify both \code{not} and \code{corner}.
##' @param center set to \code{FALSE} to not center columns of \code{X} first
##' @return list with elements \code{X, R, Rinv, xbar} where \code{xbar} is the vector of means (vector of zeros if \code{center=FALSE})
##' @examples
##' \dontrun{
##'   x <- 1 : 10
##'   X <- cbind(x, x^2)
##'   w <- selectedQr(X)
##'   w
##'   with(w, X %*% R)  # = scale(X, center=TRUE, scale=FALSE)
##'   Xqr <- w$X
##'   plot(X[, 1], Xqr[, 1])
##'   plot(X[, 1], Xqr[, 2])
##'   cov(X)
##'   cov(Xqr)
##'   X <- cbind(x, x^3, x^4, x^2)
##'   w <- selectedQr(X, not=2:3)
##'   with(w, X %*% R)
##' }
##' @author Ben Goodrich and Frank Harrell
selectedQr <- function(X, not=NULL, corner=FALSE, center=TRUE) {
  if(center) {
    X    <- scale(X, center=TRUE, scale=FALSE)
    xbar <- as.vector(attr(X, 'scaled:center'))
  } else xbar <- rep(0., ncol(X))
  
  if(length(not)) {
    if(corner) stop('may not specify both not and corner=TRUE')
    p <- ncol(X)
    # Handle the case where no variables are orthogonalized
    if(length(not) == p)
      return(list(X=X, R=diag(p), Rinv=diag(p), xbar=xbar))
    Xo <- X
    X  <- Xo[, -not, drop=FALSE]
    }
  p     <- ncol(X)
  QR    <- qr(X)
  Q     <- qr.Q(QR)
  R     <- qr.R(QR)
  sgns  <- sign(diag(R))
  X     <- sweep(Q, MARGIN = 2, STATS = sgns, FUN = `*`)
  R_ast <- sweep(R, MARGIN = 1, STATS = sgns, FUN = `*`)
  cornr <- if(corner) R_ast[p, p] else 1.
  R_ast_inverse <- backsolve(R_ast, diag(p))
  X <- X * cornr
  R_ast <- R_ast / cornr
  R_ast_inverse <- R_ast_inverse * cornr

  if(length(not)) {
    Xo[, -not]       <- X
    R <- Rinv        <- diag(p + length(not))
    R[-not,    -not] <- R_ast
    Rinv[-not, -not] <- R_ast_inverse
    return(list(X=Xo, R=R, Rinv=Rinv, xbar=xbar))
    }
  list(X = X, R = R_ast, Rinv = R_ast_inverse, xbar=xbar)
}

utils::globalVariables('est')    # why is this possibly needed?

## Code previously used with qr.stan
#  d <- list(X=X, N=N, p=p)
#  stanloc <- .Options$stancompiled
#  if(! length(stanloc)) stop('options(stancompiled) not defined')
#  m <- readRDS(paste0(stanloc, '/qr.rds'))
#  z <- rstan::sampling(m, data=d, chains=1, iter=1, algorithm='Fixed_param')
#  z <- as.matrix(z)[,]
#  n <- names(z)
#  xqr <- n[grep('Xqr\\[',  n)]
# ri  <- n[grep('Ri\\[',   n)]
#
#  ro <- as.integer(gsub('Xqr\\[(.*),.*',    '\\1', xqr))
#  co <- as.integer(gsub('Xqr\\[.*,(.*)\\]', '\\1', xqr))
#  Xqr <- matrix(NA, nrow=N, ncol=p)
#  Xqr[cbind(ro, co)] <- z[xqr]
#  ro <- as.integer(gsub('Ri\\[(.*),.*',    '\\1', ri))
#  co <- as.integer(gsub('Ri\\[.*,(.*)\\]', '\\1', ri))
#  R <- matrix(NA, nrow=p, ncol=p)
#  R[cbind(ro, co)] <- z[ri]
#  list(Xqr=Xqr, R=R)

