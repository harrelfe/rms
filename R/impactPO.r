#' Impact of Proportional Odds Assumpton
#'
#' Checks the impact of the proportional odds assumption by comparing predicted cell probabilities from a PO model with those from a multinomial or partial proportional odds logistic model that relax assumptions.  For a given model formula, fits the model with both `lrm` and either `nnet::multinom` or `VGAM::vglm` or both, and obtains predicted cell probabilities for the PO and relaxed models on the `newdata` data frame.  A `print` method formats the output.
#'
#' @param formula a model formula.  To work properly with `multinom` or `vglm` the terms should have completely specified knot locations if a spline function is being used.
#' @param relax defaults to `"both"` if `nonpo` is given, resulting in fitting two relaxed models.  Set `relax` to `"multinomial"` or `"ppo"` to fit only one relaxed model.  The multinomial model does not assume PO for any predictor.
#' @param nonpo a formula with no left hand side variable, specifying the variable or variables for which PO is not assumed.  Specifying `nonpo` results in a relaxed fit that is a partial PO model fitted with `VGAM::vglm`.  
#' @param newdata a data frame or data table with one row per covariate setting for which predictions are to be made
#' @param ... other parameters to pass to `lrm` and `multinom` such as `data=`.
#' @return an `impactPO` object which is a list with elements `estimates`, `stats`, and `mad`.  `estimates` is a data frame containing the variables and values in `newdata` in a tall and thin format with additional variable `method` ("PO", "Multinomial", "PPO"), `y` (current level of the dependent variable), and `Probability` (predicted cell probability for covariate values and value of `y` in the current row).  `stats` is a data frame containing `Deviance` the model deviance, `d.f.` the total number of parameters counting intercepts, `AIC`, `p` the number of regression coefficients, `LR chi^2` the likelihood ratio chi-square statistic for testing the predictors, `LR - p` a chance-corrected LR chi-square, `LR chi^2 test for PO` the likelihood ratio chi-square test statistic for testing the PO assumption (by comparing -2 log likelihood for a relaxed model to that of a fully PO model), `  d.f.` the degrees of freedom for this test, `  Pr(>chi^2)` the P-value for this test, `Cox-Snell R2`, `Cox-Snell R2 adj` (adjusted version of Cox-Snell R2 that is very similar to the way adjusted R2 is computed in the linear model, resulting in the regression d.f. being subtracted from the likelihood ratio chi-square statistic), `McFadden R2`, `McFadden R2 adj` (an AIC-like adjustment proposed by McFadden without full justification), `Mean |difference} from PO` the overall mean absolute difference between predicted probabilities over all categories of Y and over all covariate settings.  `mad` contains `newdata` and separately by rows in `newdata` the mean absolute difference (over Y categories) between estimated probabilities by the indicated relaxed model and those from the PO model. 
#'
#' @author Frank Harrell <fh@fharrell.com>
#' @export
#'
#' @import rms
#' @export
#' @seealso [nnet::multinom()], [VGAM::vglm()], [lrm()], [Hmisc::propsPO()]
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
                     nonpo, newdata, ...) {

  if(relax != 'multinomial' && missing(nonpo))
    stop('must specify nonpo when relax is not "multinomial"')
  if(relax != 'ppo') if(! requireNamespace('nnet', quietly=TRUE))
                       stop("This function requires the 'nnet' package")
  if(relax != 'multinomial') if(! requireNamespace('VGAM', quietly=TRUE))
                       stop("This function requires the 'VGAM' package")
  
  f <- lrm(formula, ...)
  a <- predict(f, newdata, type='fitted.ind')
  nam <- names(f$freq)
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
    r2    <- 1. - exp(- LR / n)
    r2adj <- 1. - exp(- max(0., LR - p) / n)
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
               `Cox-Snell R2`     = r2,
               `Cox-Snell R2 adj` = r2adj,
               `McFadden R2`      = 1. - dev / dev0,
               `McFadden R2 adj`  = 1. - (dev + 2 * p) / dev0,
               `Mean |difference| from PO` = if(mpo) NA else mean(abs(probs - probsPO)),
               check.names=FALSE)
    }
  stats <- st('PO', f, probsPO)

  mad <- NULL
  
  if(relax != 'multinomial') {
    ppo  <- formula(paste('FALSE ~', as.character(nonpo)[-1]))
    g <- VGAM::vglm(formula, VGAM::cumulative(parallel=ppo, reverse=TRUE), ...)
    b <- VGAM::predict(g, newdata, type='response')
    if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
    colnames(b) <- nam
    md  <- apply(abs(b - probsPO), 1, mean)
    mad <- rbind(mad, cbind(method='PPO', newdata, `Mean |difference|`=md))
    A <- rbind(A, cbind(method='PPO', newdata, b))
    stats <- rbind(stats, st('PPO', g, b))
    }

  if(relax != 'ppo') {
    g <- nnet::multinom(formula, ..., trace=FALSE)
    b <- predict(g, newdata, 'probs')
    if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
    colnames(b) <- nam
    md  <- apply(abs(b - probsPO), 1, mean)
    mad <- rbind(mad, cbind(method='Multinomial', newdata, `Mean |difference|`=md))
     A <- rbind(A, cbind(method='Multinomial', newdata, b))
    stats <- rbind(stats, st('Multinomial', g, b))
  }

  z <- reshape(A, direction='long', varying=list(nam),
               times=nam, v.names='Probability', timevar='y')
  z$method <- factor(z$method,
   c('PO', c('PPO', 'Multinomial')[c(relax != 'multinomial', relax != 'PPO')]))

  rownames(stats) <- NULL
                                    
  structure(list(estimates=z, stats=stats, mad=mad),
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

  invisible()
}
