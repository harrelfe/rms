#' Examine proportional odds and parallelism assumptions of `orm`, `lrm` and `blrm` model fits.
#'
#' Based on codes and strategies from Frank Harrell's canonical `Regression Modeling Strategies` text
#'
#' Strategy 1: Apply different link functions to Prob of Binary Ys (defined by cutval). Regress transformed outcome on combined X and assess constancy of slopes (betas) across cut-points \cr
#' Strategy 2: Generate score residual plot for each predictor (for response variable with <10 unique levels) \cr
#' Strategy 3: Assess parallelism of link function transformed inverse CDFs curves for different combined (XBeta) or specific predictor levels (for response variables with >=10 unique levels)
#'
#' @param mod.orm Model fit of class `orm`, `lrm`, or `blrm`. `poma` will fit an `orm` object from an `blrm` object. For `fit.mult.impute` objects, `poma` will refit model on a singly-imputed data-set
#' @predvar An optional character vector for which link-transformed inverse ECDF curves are generated   
#' @param cutval Numeric vector; sequence of observed values to categorize outcome
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
#'poma(mod.orm, predvar = "hp")  ## inverse ECDF curves for "hp" variable 
#'
#'
#' ## orm model using imputation
#' dat <- mtcars
#' ## introduce NAs
#' dat[sample(rownames(dat), 10), "cyl"] <- NA
#' im <- aregImpute(~ cyl + wt + mpg + am, data = dat)
#' aa <- fit.mult.impute(mpg ~ cyl + wt , xtrans = im, data = dat, fitter = orm)
#' poma(aa)



poma <- function(mod.orm, predvar=NULL,cutval) {
  
  ### Ensure that lrm and orm objects are supplied
  if(!any(class(mod.orm) %in% Cs(lrm, orm, blrm))) {
    stop('rms object must be of class (b)lrm or orm', call. = FALSE)
  }
  
  ## Get data and response variable
  ## if NAs in dataset, issue warnings and remove NAs 
  data = eval(mod.orm$call$data)[all.vars(mod.orm$sformula)]
  data = data.frame(data)  
  if (anyNA(data)) {
    warning("NA responses in sample")              ## issue warnings about NAs
    data <- data[complete.cases(data), ]  ## remove NAs (if any)
  }
  mydv <-   eval (data) [ , all.vars(mod.orm$sformula)[1] ]
  
  ## Exclude binary logistic models
  if(length(unique(mydv)) == 2L){
    stop("'poma` is not applicable to binary logistic models", call. = FALSE)
  }
  
  ## ensure valid input values for `predvar`
  predvar_valid = mod.orm$Design$name  
  if(!is.null(predvar)) {
    if(length(unique(mydv)) <= 10) {
      warning("`predvar` ignored as number of response levels <=10")  
    } else 
      if (predvar %nin% predvar_valid) {
        stop(predvar, " is not found in the model. ", "Available predictors are:\n",  paste(predvar_valid, collapse=', '))
      }
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
  
  ## Fit a (Frequentist) `orm` from a Bayesian `blrm` 
  if(any(class(mod.orm) %in% Cs(blrm))) 
    if(length(mod.orm$cppo))
      stop('poma not yet implemented for cppo models') else {
        cat("Fitting an `orm` model from a `blrm` model\n\n")
        mod.orm <-  orm(formula(mod.orm), x=TRUE, y=TRUE, data=eval(data))
      }
  
  
  ### Convert DV into numeric vector when factor DV is supplied
  cat("Unique values of Y:\n", unique(sort(mydv)), "\n")
  
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
      co <- coef(suppressWarnings(glm(mod.ols$y < k ~ combined_x, data=eval(data), family=binomial (link))))
      r <-  rbind(r, data.frame (link=link, cut.value=k, slope =round(co[2],3)))
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
    
  } else 
    if(!is.null(predvar)){
      ## inverse ECDF curves for different levels of user-specified predictors 
      
      poma_theme <-
        theme_light(base_size = 13) +
        theme(
          legend.position= c(0.7,0.7), 
          legend.text=element_text(size=10), 
          legend.key.size=unit(3, 'mm'),
          axis.title.y = element_blank(),    
          text = element_text(colour = d.gray),
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(face = "italic"),
          panel.grid.major = element_line(colour = "#E7E9E9"),
          panel.grid.minor = element_line(colour = "#F4F5F5"),
          strip.text = element_text(face = "bold", colour = "#2E3F4F"  ),
          strip.background = element_rect(colour = "#2E3F4F", fill = "grey90")      
        )
      
      isdiscrete <- function(z) is.factor(z) || is.character(z) || 
        length(unique(z[!is.na(z)])) <= 3
      
      ecdf_df <- data.frame(y = mydv, group_vector = data[,predvar])
      ecdf_df <- ecdf_df[order(ecdf_df$group),]
      ecdf_df$group_vector <- 
        if(!isdiscrete(data[,predvar]))
          cut2(ecdf_df$group_vector, g=2)  else 
            as.character(ecdf_df$group_vector)
      
      ecdf_df_tmp <- ggplot(ecdf_df, aes(x = y, color=group_vector)) + stat_ecdf() 
      ecdf_df1 <- layer_data (ecdf_df_tmp)
      ecdf_df1$logity = plogis(1-ecdf_df1$y)
      ecdf_df1$qnormy = qnorm(1-ecdf_df1$y)
      ecdf_df1$group_vector = factor(ecdf_df1$group, 
                                     labels =if(is.factor(ecdf_df$group_vector))  
                                       levels(ecdf_df$group_vector) else unique(ecdf_df$group_vector))
      
      logit_plot <- ggplot(ecdf_df1, aes(x=x, y = logity, color=group_vector)) + 
        geom_step() +
        labs(x = all.vars(mod.orm$sformula)[1]) + 
        scale_colour_manual(name=predvar, values= unique(ecdf_df1$colour)) + 
        facet_grid( . ~ "logit~(P(Y>=y~symbol('\\275')~x))", labeller=label_parsed) +
        poma_theme 
      
      qnorm_plot <- ggplot(ecdf_df1, aes(x=x, y = qnormy, color=group_vector)) + 
        geom_step() +
        labs(x = all.vars(mod.orm$sformula)[1]) + 
        scale_colour_manual(name=predvar, values= unique(ecdf_df1$colour)) + 
        facet_grid( . ~ "phi^-1~P(Y>=y~symbol('\\275')~x)", labeller=label_parsed) +
        poma_theme 
      
      allplots <- 
        arrGrob(logit_plot , qnorm_plot,layout_matrix = rbind(c(1,2),c(1,2)  ))
      
      
      return(allplots)
      
      
    } else {
      ## inverse ECDF curves for different levels of combined X values
      ## Assess parallelism of link function transformed inverse CDFs curves
      ## Codes to generate curves are from Harrell's rms book p368-369
      p <- function (fun, row, col) {
        f <-  substitute (fun)
        g <-  function (F) eval(f)
        
        ## Number of groups (2 to 5) based on sample size
        ecdfgroups = pmax(2, pmin(5, round( nrow(data)/20)))
        
        
        
        z <-  Ecdf (~ mydv,
                    groups = 
                      if(length(unique(round(combined_x))) == 2L) 
                        cut2(combined_x, mean(combined_x)) else cut2(combined_x, g = ecdfgroups),
                    fun = function (F) g(1 - F),
                    xlab = all.vars(mod.ols$sformula)[[1]],
                    ylab = as.expression (f) ,
                    xlim = c(quantile(mod.orm$y, 0.02, na.rm= TRUE), quantile(mod.orm$y, 0.98, na.rm= TRUE)),
                    label.curve= FALSE)
        print (z, split =c(col , row , 2, 2) , more = row < 2 | col < 2)
      }
      
      p (fun = log (F/(1-F)), 1, 1)
      p (fun = qnorm(F), 1, 2)
      p (fun = log (-log (1-F)), 2, 1)
      p( fun = -log (-log (F)), 2, 2)
    }
  
  
}