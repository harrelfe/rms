#' Impact of Proportional Odds Assumpton

#' Checks the impact of the proportional odds assumption by comparing predicted cell probabilities from a PO model with those from a multinomial or partial proportional odds logistic model that relax assumptions.  For a given model formula, fits the model with both `lrm` and either `nnet::multinom` or `VGAM::vglm` or both, and obtains predicted cell probabilities for the PO and relaxed models on the `newdata` data frame.
#'
#' @param formula a model formula.  To work properly with `multinom` or `vglm` the terms should have completely specified knot locations if a spline function is being used.
#' @param relax defaults to `"both"` if `nonpo` is given, resulting in fitting two relaxed models.  Set `relax` to `"multinomial"` or `"ppo"` to fit only one relaxed model.  The multinomial model does not assume PO for any predictor.
#' @param nonpo a formula with no left hand side variable, specifying the variable or variables for which PO is not assumed.  Specifying `nonpo` results in a relaxed fit that is a partial PO model fitted with `VGAM::vglm`.  
#' @param newdata a data frame or data table with one row per covariate setting for which predictions are to be made
#' @param ... other parameters to pass to `lrm` and `multinom` such as `data=`.
#' @return an `impactPO` object which is a list with elements `estimates` and `stats`.  `data` is a data frame containing the variables and values in `newdata` in a tall and thin format with additional variable `method` ("PO", "multinomial", "PPO"), `y` (current level of the dependent variable), and `Probability` (predicted cell probability for covariate values and value of `y` in the current row).
#'
#' @author Frank Harrell <fh@fharrell.com>
#'
#' @import rms
#'#' @export
#'#' @seealso [nnet::mulinom()], [VGAM::vglm()], [lrm()], [Hmisc::propsPO()]
#'#' @keywords category models regression
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
#' w   <- impactPO(y ~ age + sex, newdata=d)
#' w
#' # Note that PO model is a better model than multinomial (lower AIC)
#' # since multinomial model's improvement in fit is low in comparison
#' # with number of additional parameters estimated
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
#' d <- expand.grid(sex=c('female', 'male'), age=c(40, 60))
#' w <- impactPO(y ~ age + sex, newdata=d)
#' w
#' ggplot(w$estimates, aes(x=method, y=Probability, fill=revo(y))) +
#'   facet_grid(age ~ sex) + geom_col() +
#'  xlab('') + guides(fill=guide_legend(title=''))
#'
#' Relax PO assumption only for sex and also show multinomial estimates
#' w <- impactPO(y ~ age + sex, nonpo = ~ sex, newdata=d)
#' # use relax='ppo' to omit the multinomial model 
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
  A <- cbind(method='PO', newdata, a)
  stats <- data.frame(method   = 'PO',
                      Deviance = deviance(f)[2],
                      `d.f.`   = length(coef(f)),
                      AIC      = AIC(f))

  if(relax != 'ppo') {
    g <- nnet::multinom(formula, ..., trace=FALSE)
    b <- predict(g, newdata, 'probs')
    if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
    colnames(b) <- nam
    A <- rbind(A, cbind(method='multinomial', newdata, b))
    stats <- rbind(stats,
                   data.frame(method   = 'multinomial',
                              Deviance = deviance(g),
                              `d.f.`   = length(coef(g)),
                              AIC      = deviance(g) + 2 * length(coef(g))))
  }

  if(relax != 'multinomial') {
    ppo  <- formula(paste('FALSE ~', as.character(nonpo)[-1]))
    g <- VGAM::vglm(formula, VGAM::cumulative(parallel=ppo, reverse=TRUE), ...)
    b <- VGAM::predict(g, newdata, type='response')
    if(nrow(newdata) == 1) b <- matrix(b, nrow=1)
    colnames(b) <- nam
    A <- rbind(A, cbind(method='PPO', newdata, b))
    stats <- rbind(stats,
                   data.frame(method   = 'PPO',
                              Deviance = deviance(g),
                              `d.f.`   = length(coef(g)),
                              AIC      = deviance(g) + 2 * length(coef(g))))
    }
  
  z <- reshape(A, direction='long', varying=list(nam),
               times=nam, v.names='Probability', timevar='y')
  z$method <- factor(z$method,
   c('PO', c('PPO', 'multinomial')[c(relax != 'multinomial', relax != 'PPO')]))

  rownames(stats) <- NULL
  list(estimates=z, stats=stats)
}
