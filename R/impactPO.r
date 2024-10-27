#' Impact of Proportional Odds Assumpton
#'
#' Checks the impact of the proportional odds assumption by comparing predicted cell probabilities from a PO model with those from a multinomial or partial proportional odds logistic model that relax assumptions.  For a given model formula, fits the model with both `lrm` and either `nnet::multinom` or `VGAM::vglm` or both, and obtains predicted cell probabilities for the PO and relaxed models on the `newdata` data frame.  A `print` method formats the output.
#'
#' Since partial proportional odds models and especially multinomial logistic models can have many parameters, it is not feasible to use this model comparison approach when the number of levels of the dependent variable Y is large.  By default, the function will use [Hmisc::combine.levels()] to combine consecutive levels if the lowest frequency category of Y has fewer than `minfreq` observations.
#'
#' @param formula a model formula.  To work properly with `multinom` or `vglm` the terms should have completely specified knot locations if a spline function is being used.
#' @param relax defaults to `"both"` if `nonpo` is given, resulting in fitting two relaxed models.  Set `relax` to `"multinomial"` or `"ppo"` to fit only one relaxed model.  The multinomial model does not assume PO for any predictor.
#' @param nonpo a formula with no left hand side variable, specifying the variable or variables for which PO is not assumed.  Specifying `nonpo` results in a relaxed fit that is a partial PO model fitted with `VGAM::vglm`.  
#' @param newdata a data frame or data table with one row per covariate setting for which predictions are to be made
#' @param data data frame containing variables to fit; default is the frame in which `formula` is found
#' @param minfreq minimum sample size to allow for the least frequent category of the dependent variable.  If the observed minimum frequency is less than this, the [Hmisc::combine.levels()] function will be called to combine enough consecutive levels so that this minimum frequency is achieved.
#' @param B number of bootstrap resamples to do to get confidence intervals for differences in predicted probabilities for relaxed methods vs. PO model fits.  Default is not to run the bootstrap.  When running the bootstrap make sure that all model variables are explicitly in `data=` so that selection of random subsets of data will call along the correct rows for all predictors.
#' @param ... other parameters to pass to `lrm` and `multinom`
#' @return an `impactPO` object which is a list with elements `estimates`, `stats`, `mad`, `newdata`, `nboot`, and `boot`.  `estimates` is a data frame containing the variables and values in `newdata` in a tall and thin format with additional variable `method` ("PO", "Multinomial", "PPO"), `y` (current level of the dependent variable), and `Probability` (predicted cell probability for covariate values and value of `y` in the current row).  `stats` is a data frame containing `Deviance` the model deviance, `d.f.` the total number of parameters counting intercepts, `AIC`, `p` the number of regression coefficients, `LR chi^2` the likelihood ratio chi-square statistic for testing the predictors, `LR - p` a chance-corrected LR chi-square, `LR chi^2 test for PO` the likelihood ratio chi-square test statistic for testing the PO assumption (by comparing -2 log likelihood for a relaxed model to that of a fully PO model), `  d.f.` the degrees of freedom for this test, `  Pr(>chi^2)` the P-value for this test, `MCS R2` the Maddala-Cox-Snell R2 using the actual sample size, `MCS R2 adj` (`MCS R2` adjusted for estimating `p` regression coefficients by subtracting `p` from `LR`), `McFadden R2`, `McFadden R2 adj` (an AIC-like adjustment proposed by McFadden without full justification), `Mean |difference} from PO` the overall mean absolute difference between predicted probabilities over all categories of Y and over all covariate settings.  `mad` contains `newdata` and separately by rows in `newdata` the mean absolute difference (over Y categories) between estimated probabilities by the indicated relaxed model and those from the PO model.  `nboot` is the number of successful bootstrap repetitions, and `boot` is a 4-way array with dimensions represented by the `nboot` resamples, the number of rows in `newdata`, the number of outcome levels, and elements for `PPO` and `multinomial`.  For the modifications of the Maddala-Cox-Snell indexes see `Hmisc::R2Measures`.
#'
#' @author Frank Harrell <fh@fharrell.com>
#' @export
#'
#' @import rms
#' @export
#' @seealso [nnet::multinom()], [VGAM::vglm()], [lrm()], [Hmisc::propsPO()], [Hmisc::R2Measures()], [Hmisc::combine.levels()]
#' @keywords category models regression
#' @references
#' [Adjusted R-square note](https://hbiostat.org/bib/r2.html)
#' @md
#'
#' @examples
#'
#' \dontrun{
#' set.seed(1)
#' age <- rnorm(500, 50, 10)
#' sex <- sample(c('female', 'male'), 500, TRUE)
#' y   <- sample(0:4, 500, TRUE)
#' d   <- expand.grid(age=50, sex=c('female', 'male'))
#' w   <- impactPO(y ~ age + sex, nonpo = ~ sex, newdata=d)
#' w
#' # Note that PO model is a better model than multinomial (lower AIC)
#' # since multinomial model's improvement in fit is low in comparison
#' # with number of additional parameters estimated.  Same for PO model
#' # in comparison with partial PO model.
#'
#' # Reverse levels of y so stacked bars have higher y located higher
#' revo <- function(z) {
#'   z <- as.factor(z)
#'   factor(z, levels=rev(levels(as.factor(z))))
#' }
#'
#' require(ggplot2)
#' ggplot(w$estimates, aes(x=method, y=Probability, fill=revo(y))) +
#'   facet_wrap(~ sex) + geom_col() +
#'   xlab('') + guides(fill=guide_legend(title=''))
#'
#' # Now vary 2 predictors
#' 
#' d <- expand.grid(sex=c('female', 'male'), age=c(40, 60))
#' w <- impactPO(y ~ age + sex, nonpo = ~ sex, newdata=d)
#' w
#' ggplot(w$estimates, aes(x=method, y=Probability, fill=revo(y))) +
#'   facet_grid(age ~ sex) + geom_col() +
#'  xlab('') + guides(fill=guide_legend(title=''))
#' }

impactPO <- function(formula,
                     relax=if(missing(nonpo)) 'multinomial' else 'both',
                     nonpo, newdata, data=environment(formula),
                     minfreq=15, B=0, ...) {

  do.mn  <- relax != 'ppo'
  do.ppo <- relax != 'multinomial'
  if(do.ppo && missing(nonpo))
    stop('must specify nonpo when relax is not "multinomial"')
  if(do.mn) if(! requireNamespace('nnet', quietly=TRUE))
                       stop("This function requires the 'nnet' package")
  if(do.ppo) if(! requireNamespace('VGAM', quietly=TRUE))
                       stop("This function requires the 'VGAM' package")

  yvarname <- all.vars(formula)[1]
  yv       <- data[[yvarname]]
  freq     <- table(yv)
  if(min(freq) < minfreq) {
    message(paste('\nimpactPO: Minimum frequency of a distinct Y value is',
                  min(freq), 'which is below', minfreq,
                  'so combine.levels\nis used to pool consecutive levels until the minimum frequency exceeds',
                  minfreq, '\n'))
    yv  <- combine.levels(yv, m=minfreq, ord=TRUE)
    if(length(levels(yv)) < 3)
      stop(paste('could not get at least 3 Y levels with >=', minfreq,
                 'observations.  Lower minfreq if you want to take a risk.'))
    message(paste('New Y levels:', paste(levels(yv), collapse='; '), '\n'))
    data$.Yreduced. <- yv
    ## keeps original y variable in calling environment from being modified
    formula <- update(formula, .Yreduced. ~ .)
    }
  f <- lrm(formula, data=data, ...)
  a <- predict(f, newdata, type='fitted.ind')
  ytable <- f$freq
  nam <- names(ytable)
  if(nrow(newdata) == 1) a <- matrix(a, nrow=1)
  colnames(a) <- nam
  probsPO     <- a
  A <- cbind(method='PO', newdata, a)
  dev   <- deviance(f)
  dev0  <- dev[1]
  devPO <- dev[2]
  dfPO  <- f$stats['d.f.']
  n     <- f$stats['Obs']
  kint  <- length(nam) - 1
  st <- function(method, fit, probs) {
    df  <- length(coef(fit))
    p   <- df - kint
    dev <- deviance(fit)
    dev <- dev[length(dev)]
    aic <- dev + 2 * df
    LR  <- dev0 - dev
    r2m <- R2Measures(LR, p, n, ytable)
    r2    <- r2m[1]
    r2adj <- r2m[2]
    mpo   <- method == 'PO'
    data.frame(method     = method,
               Deviance   = dev,
               `d.f.`     = df,
               AIC        = aic,
               p          = p,
               `LR chi^2` = LR,
               `LR - p`   = LR - p,
               `LR chi^2 test for PO` = if(mpo) NA else devPO - dev,
               `  d.f.`    = if(mpo) NA else p - dfPO,
               `  Pr(>chi^2)` = if(mpo) NA else 1. - pchisq(devPO - dev, p - dfPO),
               `MCS R2`     = r2,
               `MCS R2 adj` = r2adj,
               `McFadden R2`      = 1. - dev / dev0,
               `McFadden R2 adj`  = 1. - (dev + 2 * p) / dev0,
               `Mean |difference| from PO` = if(mpo) NA else mean(abs(probs - probsPO)),
               check.names=FALSE)
    }
  stats <- st('PO', f, probsPO)

  mad <- NULL
  
  if(do.ppo) {
    ppo  <- formula(paste('FALSE ~', as.character(nonpo)[-1]))
    g <- VGAM::vglm(formula, VGAM::cumulative(parallel=ppo, reverse=TRUE),
                    data=data, ...)
    vglmcoef <- coef(g)   # save to jump start bootstrap estimates
    b <- VGAM::predict(g, newdata, type='response')
    if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
    colnames(b) <- nam
    probsPPO <- b
    md  <- apply(abs(b - probsPO), 1, mean)
    mad <- rbind(mad, cbind(method='PPO', newdata, `Mean |difference|`=md))
    A <- rbind(A, cbind(method='PPO', newdata, b))
    stats <- rbind(stats, st('PPO', g, b))
    }

  if(do.mn) {
    g <- nnet::multinom(formula, data=data, ..., trace=FALSE)
    b <- predict(g, newdata, 'probs')
    if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
    colnames(b) <- nam
    probsM      <- b
    md  <- apply(abs(b - probsPO), 1, mean)
    mad <- rbind(mad, cbind(method='Multinomial', newdata, `Mean |difference|`=md))
    A <- rbind(A, cbind(method='Multinomial', newdata, b))
    stats <- rbind(stats, st('Multinomial', g, b))
  }

  z <- reshape(A, direction='long', varying=list(nam),
               times=nam, v.names='Probability', timevar='y')
  z$method <- factor(z$method,
   c('PO', c('PPO', 'Multinomial')[c(do.ppo, do.mn)]))

  rownames(stats) <- NULL

  nboot <- 0
  boot  <- array(NA, c(B, nrow(newdata), length(nam), do.ppo + do.mn),
                 dimnames=list(NULL, NULL, nam,
                               c('PPO', 'Multinomial')[c(do.ppo, do.mn)]))
  if(B > 0) {
    if(! is.data.frame(data)) data <- model.frame(formula, data=data)
    for(i in 1 : B) {
      j   <- sample(nrow(data), nrow(data), replace=TRUE)
      dat <- data[j, ]

      f <- lrm(formula, data=dat, ...)
      if(length(f$fail) && f$fail) next
      # If a Y category was not selected in this bootstrap sample,
      # go to the next sample
      if(length(names(f$freq)) != length(nam)) next
      a <- predict(f, newdata, type='fitted.ind')
      if(nrow(newdata) == 1) a <- matrix(a, nrow=1)
      colnames(a) <- nam

 if(do.ppo) {
        g <- try(VGAM::vglm(formula,
                            VGAM::cumulative(parallel=ppo, reverse=TRUE),
                            data=dat, coefstart=vglmcoef, ...),
                 silent=TRUE)
        if(inherits(g, 'try-error')) next
        b <- VGAM::predict(g, newdata, type='response')
        if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
        colnames(b) <- nam
        pppo <- b
        }

      if(do.mn) {
        g <- try(nnet::multinom(formula, data=dat, ..., trace=FALSE))
        if(inherits(g, 'try-error')) next
        b <- predict(g, newdata, 'probs')
        if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
        colnames(b) <- nam
        pmn <- b
      }
      nboot <- nboot + 1
      if(do.ppo) boot[nboot, , , 'PPO']         <- a - pppo
      if(do.mn)  boot[nboot, , , 'Multinomial'] <- a - pmn
    }
    if(nboot < B) boot <- boot[1 : nboot, , , , drop=FALSE]
    }
                                    
  structure(list(estimates=z, stats=stats, mad=mad, newdata=newdata,
                 nboot=nboot, boot=if(B > 0) boot),
            class='impactPO')
}


##' Print Result from impactPO
##'
##' Prints statistical summaries and optionally predicted values computed by `impactPO`, transposing statistical summaries for easy reading
##' @param x an object created by `impactPO`
##' @param estimates set to `FALSE` to suppess printing estimated category probabilities.  Defaults to `TRUE` when the number of rows < 16.
##' @param ... ignored
##' @author Frank Harrell
##' @method print impactPO
##' @export
##' @md
print.impactPO <- function(x, estimates=nrow(x$estimates) < 16, ...) {
  stats <- x$stats
  fstats <- stats
  integercol <- c('p', 'd.f.', '  d.f.')
  r2col <- c('Mean |difference| from PO',
             names(stats)[grep('R2', names(stats))])
  z <- function(x, digits=0, pval=FALSE) {
    y <- if(pval) ifelse(x < 0.0001, '<0.0001', format(round(x, 4)))
    else format(if(digits == 0) x else round(x, digits))
    y[is.na(x)] <- ''
    y
    }
  pvn <- '  Pr(>chi^2)'
  for(j in integercol) fstats[[j]]   <- z(fstats[[j]])
  for(j in r2col)      fstats[[j]]   <- z(fstats[[j]], 3)
                       fstats[[pvn]] <- z(fstats[[pvn]], pval=TRUE)
  for(j in setdiff(names(fstats),
                   c('method', integercol, r2col, pvn)))
    fstats[[j]] <- z(fstats[[j]], 2)
  fstats <- t(fstats)
  colnames(fstats) <- fstats[1, ]
  fstats <- fstats[-1, ]
  print(fstats, quote=FALSE)

  if(estimates) {
    est <- x$estimates
    est$Probability <- round(est$Probability, 4)
    cat('\n')
    print(est)
  }

  cat('\nCovariate combination-specific mean |difference| in predicted probabilities\n\n')
  x$mad$`Mean |difference|` <- round(x$mad$`Mean |difference|`, 3)
  print(x$mad)

  if(x$nboot > 0) {
    boot <- x$boot
    cat('\nBootstrap 0.95 confidence intervals for differences in model predicted\nprobabilities based on', x$nboot, 'bootstraps\n\n')
    nd <- nrow(x$newdata)
    cl <- function(x) {
      qu <- unname(quantile(x, c(0.025, 0.975)))
      c(Lower=qu[1], Upper=qu[2]) }
    for(i in 1 : nd) {
      cat('\n')
      print(x$newdata[i, ])
      b <- boot[, i, , , drop=FALSE]
      b <- round(apply(b, 2 : 4, cl), 3)
      for(model in dimnames(b)[[4]]) {
        cat('\nPO - ', model, ' probability estimates\n\n', sep='')
        print(b[, , , model])
        }
      }
    }

  invisible()
}
