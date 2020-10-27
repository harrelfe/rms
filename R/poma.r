#' Examine proportional odds and parallelism assumptions of `orm` and `lrm` model fits.
#'
#' Based on codes and strategies from Frank Harrell's canonical `Regression Modeling Strategies` text
#'
#' Strategy 1: Apply different link functions to Prob of Binary Ys (defined by cutval). Regress transformed outcome on combined X and assess constancy of slopes (betas) across cut-points \cr
#' Strategy 2: Generate score residual plot for each predictor (for response variable with <10 unique levels) \cr
#' Strategy 3: Assess parallelism of link function transformed inverse CDFs curves for different XBeta levels (for response variables with >=10 unique levels)
#'
#' @param mod.orm Model fit of class `orm` or `lrm`. For `fit.mult.impute` objects, `poma` will refit model on a singly-imputed data-set
#' @param cutval Numeric vector; sequence of observed values to cut outcome
#'
#' @author Yong Hao Pua <puayonghao@gmail.com>
#'
#' @import rms
#'
#' @export
#'
#' @seealso Harrell FE. *Regression Modeling Strategies: with applications to linear models,
#' logistic and ordinal regression, and survival analysis.* New York: Springer Science, LLC, 2015.
#'
#' @examples
#'
#'## orm model (response variable has fewer than 10 unique levels)
#'mod.orm <- orm(carb ~ cyl + hp , x=TRUE, y=TRUE, data = mtcars)
#'poma(mod.orm)
#'
#'
#'## orm model (response variable has >=10 unique levels)
#'mod.orm <- orm(mpg ~ cyl + hp , x=TRUE, y=TRUE, data = mtcars)
#'poma(mod.orm)
#'
#'
#' ## orm model using imputation
#' dat <- mtcars
#' ## introduce NAs
#' dat[sample(rownames(dat), 10), "cyl"] <- NA
#' im <- aregImpute(~ cyl + wt + mpg + am, data = dat)
#' aa <- fit.mult.impute(mpg ~ cyl + wt , xtrans = im, data = dat, fitter = orm)
#' poma(aa)



poma <- function(mod.orm, cutval) {

  ### Ensure that lrm and orm objects are supplied
  if(!any(class(mod.orm) %in% Cs(lrm, orm))) {
    stop('rms object must be of class lrm or orm', call. = FALSE)
  }

  ## (Re-)create mod.orm from a singly-imputed dataset
  if(any(class(mod.orm) %in% Cs(fit.mult.impute) )) {
    cat("Refitting model on a singly-imputed dataset \n")
    fit_mult_call <-  as.character(mod.orm$call)
    myformula <-      fit_mult_call[[2]]
    myfitter <-       fit_mult_call[[3]]
    myaregimpute <-   fit_mult_call[[4]]
    mydata <-         fit_mult_call[[5]]
    # extract imputed values
    imputed <- impute.transcan(x = get(myaregimpute), imputation = 1, data = get(mydata),  list.out = TRUE, pr = FALSE)
    # create one imputed dataset
    imputed_df <- get(mydata)
    imputed_df[names(imputed)] <- imputed
    # recreate model
    mod.orm <- eval(parse(text = sprintf(" %s(%s, x = T, y = T, data = imputed_df)", myfitter, myformula)))
  }


  ## Generate dataset with no missing values
  data = mod.orm$call$data
  data = eval (data) [ , all.vars(mod.orm$sformula)]
  data <- data [complete.cases(data), ]


  ### Convert DV into numeric vector when factor DV is supplied
  mydv <-   eval (data) [ , all.vars(mod.orm$sformula)[1] ]
  cat("Unique values of Y:", unique(sort(mydv)), "\n")

  ### Compute combined predictor (X) values
  if(any(class(mydv) %in% "factor") ) {
    aa <- paste0("as.numeric(", mod.orm$terms[[2]], ") ~")
    rhs <- mod.orm$terms[[3]]
    bb <- paste(deparse(rhs), collapse = "")
    newformula <- paste(aa, bb)
    cat("Formula used with non-numeric DV:", newformula, "\n")
    cat("Cut-point for factor DV refers to the jth levels - not observed Y values \n")
    mod.ols <- ols(as.formula(newformula) , x=TRUE, y=TRUE, data=eval(data))
    } else {
    cat("Cut-point for continuous DV refers to observed Y values \n")
    mod.ols <-  ols(formula(mod.orm), x=TRUE, y=TRUE, data=eval(data))
    }

  combined_x <- fitted(mod.ols)


  ### Set cutpoint values for Y
  ### for factor DV: cutpoints = 2 to max no. of levels (jth levels)
  ### for continuous DV: cutpoints = y unique values (quartiles for truly continuous response var)

  if (missing(cutval)) {
    if (any(class(mydv) %in% "factor")) cutval <- seq(2, length(unique(mydv)))
    else if( length(unique(mydv)) <= 10 ) cutval <- unique(sort(mydv))[-1]
    else cutval <- quantile(unique(mydv), c(0.25, 0.5, 0.75), na.rm = TRUE)  ## quartiles as cutpoints for continuous DV
  }


  ### Apply link functions to Prob of Binary Y (defined by cutval)
  ### Regress transformed outcome as a function of combined X. Check for constancy of slopes
  ### Codes taken from rms p368
  r <- NULL
  for (link in c("logit", "probit", "cloglog")) {
    for (k in cutval) {
      co <- coef(glm(mod.ols$y < k ~ combined_x, data=eval(data), family=binomial (link)))
      r <-  rbind(r, data.frame (link=link, cut.value=k, slope =round(co[2],2)))
    }
  }
  cat("rms-368: glm cloglog on Prob[Y<j] is log-log on Prob{Y>=j} \n")
  print(r, row.names=FALSE)


  ### Graphical Assessment
  if(length(unique(mod.orm$y)) < 10) {
    par(ask=TRUE)
    ## Generate Score residual plot for each predictor/terms
    ## Adjust par(mfrow) settings based on number of terms (codes are a little unwieldy)
    numpred <- dim(mod.orm$x)[[2]]
    if(numpred >= 9 )    par(mfrow = c(3,3))
    else if (numpred >= 5) par(mfrow = c(2,3))
    else if (numpred >= 3) par(mfrow = c(2,2))
    else if (numpred >=2) par(mfrow = c(1,2))
    else par(mfrow = c(1,1))
    resid(mod.orm, "score.binary", pl=TRUE)
    par(ask=F)
    par(mfrow = c(1,1))

  } else {

    ## Assess parallelism of link function transformed inverse CDFs curves
    ## Codes to generate curves are from Harrell's rms book p368-369
    p <- function (fun, row, col) {
      f <-  substitute (fun)
      g <-  function (F) eval(f)

      ## Number of groups (2 to 5) based on sample size
      ecdfgroups = pmax(2, pmin(5, round( dim(mod.orm$x)[[1]]/20)))

      z <-  Ecdf (~ mod.ols$y,
                  groups = cut2(combined_x, g = ecdfgroups),
                  fun = function (F) g(1 - F),
                  xlab = all.vars(mod.ols$sformula)[[1]],
                  ylab = as.expression (f) ,
                  xlim = c(quantile(mod.orm$y, 0.10, na.rm= TRUE), quantile(mod.orm$y, 0.85, na.rm= TRUE)),
                  label.curve= FALSE)
      print (z, split =c(col , row , 2, 1) , more = row < 2 | col < 2)
    }

    p (fun = log (F/(1-F)), 1, 1)
    p( fun = -log ( -log (F)), 1, 2)
  }


}





