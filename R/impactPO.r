#' Impact of Proportional Odds Assumpton

#' Checks the impact of the proportional odds assumption by comparing predicted cell probabilities from a PO model with those from a multinomial or partial proportional odds logistic model.  For a given model formula, fits the model with both `lrm` and either `nnet::multinom` or `VGAM::vglm` and obtains predicted cell probabilities for both on the `newdata` data frame.
#'
#' @param formula a model formula.  To work properly with `multinom` or `vglm` the terms should have completely specified knot locations if a spline function is being used.
#' @param nonppo a formula with no left hand side variable, specifying the variable or variables for which PO is not assumed.  Specifying `nonppo` results in the alternative (relaxed) fit being a partial PO model fitted with `VGAM::vglm`.  If `nonppo` is omitted, the multimonial model is fitted and PO is not assumed for any variables.
#' @param newdata a data frame or data table with one row per covariate setting for which predictions are to be made
#' @param ... other parameters to pass to `lrm` and `multinom` such as `data=`.
#' @return an `impactPO` object which is a list with elements `estimates` and `stats`.  `data` is a data frame containing the variables and values in `newdata` in a tall and thin format with additional variable `method` ("PO", "multinomial", or "PPO"), `y` (current level of the dependent variable), and `Probability` (predicted cell probability for covariate values and value of `y` in the current row).
#'
#' @author Frank Harrell <fh@fharrell.com>
#'
#' @import rms
#'#' @export
#'#' @seealso [nnet::mulinom()], [lrm()], [Hmisc::propsPO()]
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
#' Relax PO assumption only for sex
#' w <- impactPO(y ~ age + sex, nonppo = ~ sex, newdata=d)
#' }

impactPO <- function(formula, nonppo, newdata, ...) {
  model <- if(missing(nonppo)) 'multinomial' else 'PPO'
  switch(model,
         multinomial = if(! requireNamespace('nnet', quietly=TRUE))
                         stop("This function requires the 'nnet' package"),
         PPO         = if(! requireNamespace('VGAM', quietly=TRUE))
                         stop("This function requires the 'VGAM' package"))
  
  f <- lrm(formula, ...)
  a <- predict(f, newdata, type='fitted.ind')

  if(model == 'multinomial') {
    g <- multinom(formula, ..., trace=FALSE)
    b <- predict(g, newdata, 'probs')
    aic <- AIC(g)
  } else {
    ppo  <- formula(paste('FALSE ~', as.character(nonppo)[-1]))
    g <- VGAM::vglm(formula, VGAM::cumulative(parallel=ppo, reverse=TRUE), ...)
    b <- VGAM::predict(g, newdata, type='response')
    aic <- deviance(g) + 2 * length(coef(g))
    }
  
  nam <- names(f$freq)
  if(nrow(newdata) == 1) {
    a <- matrix(a, nrow=1)
    b <- matrix(b, nrow=1)
    }
  colnames(a) <- colnames(b) <- nam

  A <- cbind(method='PO',  newdata, a)
  B <- cbind(method=model, newdata, b)
  w <- rbind(A, B)

  z <- reshape(w, direction='long', varying=list(nam),
               times=nam, v.names='Probability', timevar='y')

  k <- coef(g)
  r <- data.frame(method   = c('PO', model),
                  Deviance = c(deviance(f)[2], deviance(g)),
                  `d.f.`   = c(length(coef(f)), length(coef(g))),
                  AIC      = c(AIC(f), aic))
  rownames(r) <- NULL
  list(estimates=z, stats=r)
}
