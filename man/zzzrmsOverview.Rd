\name{rmsOverview}
\alias{rmsOverview}
\alias{rms.Overview}
\title{Overview of rms Package}

\description{
  rms is the package that goes along with the book Regression Modeling
  Strategies.  rms does regression modeling, testing, estimation,
  validation, graphics, prediction, and typesetting by storing enhanced
  model design attributes in the fit.  rms is a re-written version of the
  Design package that has improved graphics and duplicates very little
  code in the survival package.

  The package is a collection of about 180 functions that assist and
  streamline modeling, especially for biostatistical and epidemiologic
  applications.  It also contains functions for binary and ordinal
  logistic regression models and the Buckley-James multiple regression
  model for right-censored responses, and implements penalized maximum
  likelihood estimation for logistic and ordinary linear models.  rms
  works with almost any regression model, but it was especially written
  to work with logistic regression, Cox regression, accelerated failure
  time models, ordinary linear models, the Buckley-James model,
  generalized lease squares for longitudinal data (using the nlme
  package), generalized linear models, and quantile regression (using
  the quantreg package).  rms requires the Hmisc package to
  be installed.  Note that Hmisc has several functions useful for data
  analysis (especially data reduction and imputation).

  Older references below pertaining to the Design package are relevant to rms.
}

\section{Statistical Methods Implemented}{
  \itemize{
    \item Ordinary linear regression models
    \item Binary and ordinal logistic models (proportional odds
      and continuation ratio models)
    \item Cox model
    \item Parametric survival models in the accelerated failure
      time class
    \item Buckley-James least-squares linear regression model
    with possibly right-censored responses
	\item Generalized linear model
	\item Quantile regression
	\item Generalized least squares
    \item Bootstrap model validation to obtain unbiased
      estimates of model performance without requiring a
      separate validation sample
    
    \item Automatic Wald tests of all effects in the model that
      are not parameterization-dependent (e.g., tests of
      nonlinearity of main effects when the variable does
      not interact with other variables, tests of
      nonlinearity of interaction effects, tests for
      whether a predictor is important, either as a main
      effect or as an effect modifier)
	  
    \item Graphical depictions of model estimates (effect
      plots, odds/hazard ratio plots, nomograms that
      allow model predictions to be obtained manually even
      when there are nonlinear effects and interactions
      in the model)
    
    \item Various smoothed residual plots, including some new
      residual plots for verifying ordinal logistic model
      assumptions
    \item Composing S functions to evaluate the linear
      predictor (\eqn{X\hat{beta}}{X*beta hat}), hazard function, survival
      function, quantile functions analytically from the
      fitted model
    
      \item Typesetting of fitted model using LaTeX
	  
    \item Robust covariance matrix estimation (Huber or
      bootstrap)
    
    \item Cubic regression splines with linear tail restrictions (natural splines)
    \item Tensor splines
    \item Interactions restricted to not be doubly nonlinear
    \item Penalized maximum likelihood estimation for ordinary
      linear regression and logistic regression models.
      Different parts of the model may be penalized by
      different amounts, e.g., you may want to penalize
      interaction or nonlinear effects more than main
      effects or linear effects
    
    \item Estimation of hazard or odds ratios in presence of
      nolinearity and interaction
    
    \item Sensitivity analysis for an unmeasured binary confounder in a
      binary logistic model
    }
}

\section{Motivation}{
  rms was motivated by the following needs:
  \itemize{
    \item need to automatically print interesting Wald tests that can be
    constructed from the design
    \itemize{
      \item tests of linearity with respect to each predictor
      \item tests of linearity of interactions
      \item pooled interaction tests (e.g., all interactions involving race)
      \item pooled tests of effects with higher order effects
      \itemize{
        \item test of main effect not meaningful when effect in interaction
        \item pooled test of main effect + interaction effect is meaningful
        \item test of 2nd-order interaction + any 3rd-order interaction containing
        those factors is meaningful
      }
    }

    \item need to store transformation parameters with the fit
    \itemize{
      \item example: knot locations for spline functions
      \item these are "remembered" when getting predictions, unlike standard
      S or \R
      \item  for categorical predictors, save levels so that same dummy variables
      will be generated for predictions; check that all levels in out-of-data
      predictions were present when model was fitted
    }

    \item need for uniform re-insertion of observations deleted because of NAs
    when using \code{predict} without \code{newdata} or when using
    \code{resid}

    \item need to easily plot the regression effect of any predictor
    \itemize{
      \item example: age is represented by a linear spline with knots at 40 and 60y
      plot effect of age on log odds of disease, adjusting
      interacting factors to easily specified constants
      \item  vary 2 predictors: plot x1 on x-axis, separate curves for discrete
      x2 or 3d perspective plot for continuous x2
      \item if predictor is represented as a function in the model, plots
      should be with respect to the original variable:\cr
      \code{f <- lrm(y ~ log(cholesterol)+age)} \cr
      \code{plot(Predict(f, cholesterol))   # cholesterol on x-axis, default range}
    }

    \item need to store summary of distribution of predictors with the fit
    \itemize{
      \item plotting limits (default: 10th smallest, 10th largest values or \%-tiles)
      \item effect limits   (default: .25 and .75 quantiles for continuous vars.)
      \item  adjustment values for other predictors (default: median for continuous
      predictors, most frequent level for categorical ones)
      \item discrete numeric predictors: list of possible values
      example: x=0,1,2,3,5 -> by default don't plot prediction at x=4
      \item values are on the inner-most variable, e.g. cholesterol, not log(chol.)
      \item allows estimation/plotting long after original dataset has been deleted
      \item  for Cox models, underlying survival also stored with fit, so original
      data not needed to obtain predicted survival curves
    }

    \item need to automatically print estimates of effects in presence of non-
    linearity and interaction
    \itemize{
      \item example: age is quadratic, interacting with sex
      default effect is inter-quartile-range hazard ratio (for
      Cox model), for sex=reference level
      \item user-controlled effects: \code{summary(fit, age=c(30,50),
        sex="female")} -> odds ratios for logistic model, relative survival time
      for accelerated failure time survival models
      \item effects for all variables (e.g. odds ratios) may be plotted with
      multiple-confidence-level bars
    }

    \item need for prettier and more concise effect names in printouts,
    especially for expanded nonlinear terms and interaction terms
    \itemize{
      \item use inner-most variable name to identify predictors
      \item e.g. for \code{pmin(x^2-3,10)} refer to factor with legal S-name
      \code{x}
    }

    \item need to recognize that an intercept is not always a simple
    concept
    \itemize{  
      \item some models (e.g., Cox) have no intercept
      \item some models (e.g., ordinal logistic) have multiple intercepts
    }

    \item need for automatic high-quality printing of fitted mathematical
    model (with dummy variables defined, regression spline terms
    simplified, interactions "factored").  Focus is on regression splines
    instead of nonparametric smoothers or smoothing splines, so that
    explicit formulas for fit may be obtained for use outside S.
    rms can also compose S functions to evaluate \eqn{X\beta}{X*Beta} from
    the fitted model analytically, as well as compose SAS code to
    do this.

    \item need for automatic drawing of nomogram to represent the fitted model

    \item need for automatic bootstrap validation of a fitted model, with
    only one S command (with respect to calibration and discrimination)

    \item need for robust (Huber sandwich) estimator of covariance matrix,
    and be able to do all other analysis (e.g., plots, C.L.) using the
    adjusted covariances

    \item need for robust (bootstrap) estimator of covariance matrix, easily
    used in other analyses without change

    \item need for Huber sandwich and bootstrap covariance matrices adjusted
    for cluster sampling

    \item need for routine reporting of how many observations were deleted
    by missing values on each predictor (see \code{na.delete} in Hmisc)

    \item need for optional reporting of descriptive statistics for Y stratified
    by missing status of each X (see na.detail.response)

    \item need for pretty, annotated survival curves, using the same commands
    for parametric and Cox models

    \item need for ordinal logistic model (proportional odds model, continuation
    ratio model)

	\item need for estimating and testing general contrasts without having to
  be conscious of variable coding or parameter order
  }
}

\details{
  To make use of automatic typesetting features you must
  have LaTeX or one of its variants installed.\cr

  Some aspects of rms (e.g., \code{latex}) will not work correctly if
  \code{options(contrasts=)} other than \code{c("contr.treatment",
    "contr.poly")} are used.

  rms relies on a wealth of survival analysis
  functions written by Terry Therneau of Mayo Clinic.
  Front-ends have been written for several of
  Therneau's functions, and other functions have been
  slightly modified.
}

\section{Fitting Functions Compatible with rms}{
  rms will work with a wide variety of fitting
  functions, but it is meant especially for the
  following:
  \tabular{lll}{
    \bold{Function} \tab \bold{Purpose} \tab  \bold{Related S}\cr
    \tab                \tab  \bold{Functions}\cr
    \bold{\code{ols}}         \tab Ordinary least squares linear model     \tab \code{lm}\cr
    \bold{\code{lrm}}         \tab Binary and ordinal logistic regression  \tab \code{glm}\cr
    \tab model                                   \tab \code{cr.setup}\cr
    \bold{\code{orm}}         \tab Ordinal regression model                \tab \code{lrm}\cr
    \bold{\code{psm}}         \tab Accelerated failure time parametric     \tab \code{survreg}\cr
    \tab survival model                          \tab \cr
    \bold{\code{cph}}         \tab Cox proportional hazards regression     \tab \code{coxph}\cr
		\bold{\code{npsurv}}      \tab Nonparametric survival estimates \tab
\code{survfit.formula} \cr
    \bold{\code{bj}}          \tab Buckley-James censored least squares    \tab \code{survreg}\cr
    \tab linear model                            \tab \cr
    \bold{\code{Glm}}        \tab Version of \code{glm} for use with rms \tab \code{glm}\cr
    \bold{\code{Gls}}        \tab Version of \code{gls} for use with rms \tab \code{gls}\cr
	\bold{\code{Rq}}         \tab Version of \code{rq} for use with rms  \tab \code{rq}\cr
  }
}

\section{Methods in rms}{
  The following generic functions work with fits with rms in effect:
  \tabular{lll}{
    \bold{Function}           \tab  \bold{Purpose} \tab \bold{Related}\cr
    \tab                 \tab \bold{Functions}\cr
    \bold{\code{print}}       \tab Print parameters and statistics of fit \tab \cr
    \bold{\code{coef}}        \tab Fitted regression coefficients  \tab \cr
    \bold{\code{formula}}     \tab Formula used in the fit \tab \cr
    \bold{\code{specs}}       \tab Detailed specifications of fit \tab \cr
    \bold{\code{robcov}}      \tab Robust covariance matrix estimates \tab \cr
    \bold{\code{bootcov}}     \tab Bootstrap covariance matrix estimates \tab \cr
    \bold{\code{summary}}     \tab Summary of effects of predictors \tab \cr
    \bold{\code{plot.summary}} \tab Plot continuously shaded confidence \tab \cr
    \tab bars for results of summary  \tab \cr
    \bold{\code{anova}}       \tab Wald tests of most meaningful hypotheses \tab \cr
    \bold{\code{contrast}}    \tab General contrasts, C.L., tests           \tab \cr
    \bold{\code{plot.anova}}  \tab Depict results of anova graphically      \tab \code{dotchart}     \cr
    \bold{\code{Predict}}     \tab Partial predictor effects \tab  \code{predict} \cr
	\bold{\code{plot.Predict}}\tab Plot predictor effects using lattice graphics \tab \code{predict} \cr
	\bold{\code{bplot}}       \tab 3-D plot of effects of varying two \tab \cr
                              \tab continuous predictors \tab \code{image, persp, contour} \cr
    \bold{\code{gendata}}     \tab Generate data frame with predictor       \tab \code{expand.grid} \cr
    \tab combinations (optionally interactively) \tab \cr
    \bold{\code{predict}}     \tab Obtain predicted values or design matrix \tab \cr
    \bold{\code{fastbw}}      \tab Fast backward step-down variable            \tab \code{step} \cr
    \tab selection \tab \cr
    \bold{\code{residuals}}   \tab Residuals, influence statistics from fit \tab \cr
    (or \bold{\code{resid}})  \tab                         \tab \cr
    \bold{\code{which.influence}} 
    \tab Which observations are overly               \tab \code{residuals} \cr
    \tab influential \tab \cr
    \bold{\code{sensuc}}      \tab Sensitivity of one binary predictor in \tab \cr
    \tab lrm and cph models to an unmeasured \tab \cr
    \tab binary confounder \tab \cr
    \bold{\code{latex}}       \tab LaTeX representation of fitted              \tab \cr
    \tab model or \code{anova} or \code{summary} table \tab \cr
    \bold{\code{Function}}    \tab S function analytic representation          \tab \code{Function.transcan} \cr
    \tab of a fitted regression model (\eqn{X\beta}{X*Beta}) \tab \cr
    \bold{\code{hazard}}      \tab S function analytic representation          \tab \code{rcspline.restate} \cr
    \tab of a fitted hazard function (for \code{psm}) \tab \cr
    \bold{\code{Survival}}    \tab S function analytic representation of \tab \cr
    \tab fitted survival function (for \code{psm,cph}) \tab \cr
    \bold{\code{Quantile}}    \tab S function analytic representation of \tab \cr
    \tab fitted function for quantiles of \tab \cr
    \tab survival time (for \code{psm, cph}) \tab \cr
    \bold{\code{nomogram}}    \tab Draws a nomogram for the fitted model       \tab \code{latex, plot} \cr
    \bold{\code{survest}}     \tab Estimate survival probabilities             \tab \code{survfit} \cr
    \tab  (for \code{psm, cph}) \tab \cr
    \bold{\code{survplot}}    \tab Plot survival curves (psm, cph, npsurv)             \tab plot.survfit \cr
    \bold{\code{validate}}    \tab Validate indexes of model fit using         \tab val.prob \cr
    \tab resampling \tab \cr
    \bold{\code{calibrate}}   \tab Estimate calibration curve for model \tab \cr
    \tab using resampling \tab \cr
    \bold{\code{vif}}         \tab Variance inflation factors for a fit \tab \cr
    \bold{\code{naresid}}     \tab Bring elements corresponding to missing  \tab \cr
    \tab data back into predictions and residuals \tab \cr
    \bold{\code{naprint}}     \tab Print summary of missing values \tab \cr
    \bold{\code{pentrace}}    \tab Find optimum penality for penalized MLE \tab \cr
    \bold{\code{effective.df}}
    \tab Print effective d.f. for each type of  \tab \cr
    \tab variable in model, for penalized fit or  \tab \cr
    \tab pentrace result \tab \cr
    \bold{\code{rm.impute}}   \tab Impute repeated measures data with     \tab \code{transcan}, \cr
    \tab non-random dropout \tab \code{fit.mult.impute} \cr
    \tab \emph{experimental, non-functional} \tab
  }
}

\section{Background for Examples}{
  The following programs demonstrate how the pieces of
  the rms package work together.  A (usually)
  one-time call to the function \code{datadist} requires a
  pass at the entire data frame to store distribution
  summaries for potential predictor variables.  These
  summaries contain (by default) the .25 and .75
  quantiles of continuous variables (for estimating
  effects such as odds ratios), the 10th smallest and
  10th largest values (or .1 and .9 quantiles for small
  \eqn{n}) for plotting ranges for estimated curves, and the
  total range.  For discrete numeric variables (those
  having \eqn{\leq 10}{<=10} unique values), the list of unique values
  is also stored.  Such summaries are used by the
  \code{summary.rms, Predict}, and \code{nomogram.rms}
  functions.  You may save time and defer running
  \code{datadist}.  In that case, the distribution summary
  is not stored with the fit object, but it can be
  gathered before running \code{summary} or \code{plot}.

  \code{d <- datadist(my.data.frame) # or datadist(x1,x2)}\cr
  \code{options(datadist="d")        # omit this or use options(datadist=NULL)}\cr
  \code{                             # if not run datadist yet}\cr
  \code{cf <- ols(y ~ x1 * x2)}\cr
  \code{anova(f)}\cr
  \code{fastbw(f)}\cr
  \code{Predict(f, x2)}
  \code{predict(f, newdata)}

  In the \bold{Examples} section there are three detailed examples using a
  fitting function 
  designed to be used with rms, \code{lrm} (logistic
  regression model).  In \bold{Detailed Example 1} we
  create 3 predictor variables and a two binary response
  on 500 subjects.  For the first binary response, \code{dz},
  the true model involves only \code{sex} and \code{age}, and there is
  a nonlinear interaction between the two because the log
  odds is a truncated linear relationship in \code{age} for
  females and a quadratic function for males.  For the
  second binary outcome, \code{dz.bp}, the true population model
  also involves systolic blood pressure (\code{sys.bp}) through
  a truncated linear relationship.  First, nonparametric
  estimation of relationships is done using the Hmisc
  package's \code{plsmo} function which uses \code{lowess} with outlier
  detection turned off for binary responses.  Then
  parametric modeling is done using restricted cubic
  splines.  This modeling does not assume that we know
  the true transformations for \code{age} or \code{sys.bp} but that
  these transformations are smooth (which is not actually
  the case in the population).

  For \bold{Detailed Example 2}, suppose that a
  categorical variable treat has values \code{"a", "b"}, and
  \code{"c"}, an ordinal variable \code{num.diseases} has values
  0,1,2,3,4, and that there are two continuous variables,
  \code{age} and \code{cholesterol}.  \code{age} is fitted with a restricted
  cubic spline, while \code{cholesterol} is transformed using
  the transformation \code{log(cholesterol - 10)}.  Cholesterol
  is missing on three subjects, and we impute these using
  the overall median cholesterol.  We wish to allow for
  interaction between \code{treat} and \code{cholesterol}.  The
  following S program will fit a logistic model,
  test all effects in the design, estimate effects, and
  plot estimated transformations. The fit for
  \code{num.diseases} really considers the variable to be a
  5-level categorical variable. The only difference is
  that a 3 d.f. test of linearity is done to assess
  whether the variable can be re-modeled "asis".  Here
  we also show statements to attach the rms package
  and store predictor characteristics from datadist.

  \bold{Detailed Example 3} shows some of the survival
  analysis capabilities of rms related to the Cox
  proportional hazards model.  We simulate data for 2000
  subjects with 2 predictors, \code{age} and \code{sex}.  In the true
  population model, the log hazard function is linear in
  \code{age} and there is no \code{age} \eqn{\times}{x} \code{sex}
  interaction.  In the  
  analysis below we do not make use of the linearity in
  age.  rms makes use of many of Terry Therneau's
  survival functions that are builtin to S.

  The following is a typical sequence of steps that
  would be used with rms in conjunction with the Hmisc
  \code{transcan} function to do single imputation of all NAs in the
  predictors (multiple imputation would be better but would be
  harder to do in the context of bootstrap model validation),
  fit a model, do backward stepdown to reduce the number of
  predictors in the model (with all the severe problems this can
  entail), and use the bootstrap to validate this stepwise model,
  repeating the variable selection for each re-sample.  Here we
  take a short cut as the imputation is not repeated within the
  bootstrap.

  In what follows we (atypically) have only 3
  candidate predictors.  In practice be sure to have the
  validate and calibrate functions operate on a model fit that
  contains all predictors that were involved in previous analyses
  that used the response variable.  Here the imputation
  is necessary because backward stepdown would otherwise delete
  observations missing on any candidate variable.

  Note that you would have to define \code{x1, x2, x3, y} to run
  the following code.

  \code{xt <- transcan(~ x1 + x2 + x3, imputed=TRUE)}\cr
  \code{impute(xt)  # imputes any NAs in x1, x2, x3}\cr
  \code{# Now fit original full model on filled-in data}\cr
  \code{f <- lrm(y ~ x1 + rcs(x2,4) + x3, x=TRUE, y=TRUE) #x,y allow boot.}\cr
  \code{fastbw(f)}\cr
  \code{# derives stepdown model (using default stopping rule)}\cr
  \code{validate(f, B=100, bw=TRUE) # repeats fastbw 100 times}\cr
  \code{cal <- calibrate(f, B=100, bw=TRUE)  # also repeats fastbw}\cr
  \code{plot(cal)}
}

\examples{
## To run several comprehensive examples, run the following command
\dontrun{
demo(all, 'rms')
}
}

\section{Common Problems to Avoid}{
  \enumerate{
    \item Don't have a formula like \code{y ~ age + age^2}.
    In S you need to connect related variables using
    a function which produces a matrix, such as \code{pol} or
    \code{rcs}.
    This allows effect estimates (e.g., hazard ratios)
    to be computed as well as multiple d.f. tests of
    association.

    \item Don't use \code{poly} or \code{strata} inside formulas used in
    rms.  Use \code{pol} and \code{strat} instead.

    \item Almost never code your own dummy variables or
    interaction variables in S.  Let S do this
    automatically.  Otherwise, \code{anova} can't do its
    job.

    \item Almost never transform predictors outside of
    the model formula, as then plots of predicted
    values vs. predictor values, and other displays,
    would not be made on the original scale.  Use
    instead something like \code{y ~ log(cell.count+1)},
    which will allow \code{cell.count} to appear on
    \eqn{x}-axes.  You can get fancier, e.g.,
    \code{y ~ rcs(log(cell.count+1),4)} to fit a restricted
    cubic spline with 4 knots in \code{log(cell.count+1)}.
    For more complex transformations do something
    like %\cr
    \code{f <- function(x) \{}\cr
    \code{\ldots various 'if' statements, etc.}\cr
    \code{log(pmin(x,50000)+1)}\cr
    \code{\}}\cr
    \code{fit1 <- lrm(death ~ f(cell.count))}\cr
    \code{fit2 <- lrm(death ~ rcs(f(cell.count),4))}\cr
    \code{\}}

    \item Don't put \code{$} inside variable names used in formulas.
    Either attach data frames or use \code{data=}.

    \item Don't forget to use \code{datadist}.  Try to use it
    at the top of your program so that all model fits
    can automatically take advantage if its
    distributional summaries for the predictors.

    \item Don't \code{validate} or \code{calibrate} models which were
    reduced by dropping "insignificant" predictors.
    Proper bootstrap or cross-validation must repeat
    any variable selection steps for each re-sample.
    Therefore, \code{validate} or \code{calibrate} models
    which contain all candidate predictors, and if
    you must reduce models, specify the option
    \code{bw=TRUE} to \code{validate} or \code{calibrate}.

    \item Dropping of "insignificant" predictors ruins much
    of the usual statistical inference for
    regression models (confidence limits, standard
    errors, \eqn{P}-values, \eqn{\chi^2}{chi-squares}, ordinary indexes of
    model performance) and it also results in models
    which will have worse predictive discrimination.
  }
}

\section{Accessing the Package}{
  Use \code{require(rms)}.
}

\references{
  The primary resource for the rms package is
  \emph{Regression Modeling Strategies} by
  FE Harrell (Springer-Verlag, 2001) and the web page
  \url{http://biostat.mc.vanderbilt.edu/rms}.  See also
  the Statistics in Medicine articles by Harrell \emph{et al} listed
  below for case studies of modeling and model validation using rms.
  Also see the free book by Alzola and Harrell at
  \url{http://biostat.mc.vanderbilt.edu}.

  Several datasets useful for multivariable modeling with
  rms are found at
  \url{http://biostat.mc.vanderbilt.edu/DataSets}.
}

\section{Published Applications of rms and Regression Splines}{
  \itemize{
    \item Spline fits
    \enumerate{
      \item Spanos A, Harrell FE, Durack DT (1989): Differential
      diagnosis of acute meningitis: An analysis of the
      predictive value of initial observations.  \emph{JAMA}
      2700-2707.

      \item Ohman EM, Armstrong PW, Christenson RH, \emph{et al}. (1996):
      Cardiac troponin T levels for risk stratification in
      acute myocardial ischemia.  \emph{New Eng J Med} 335:1333-1341.
    }

    \item Bootstrap calibration curve for a parametric survival
    model:
    \enumerate{
      \item Knaus WA, Harrell FE, Fisher CJ, Wagner DP, \emph{et al}.
      (1993):  The clinical evaluation of new drugs for
      sepsis: A prospective study design based on survival
      analysis.  \emph{JAMA} 270:1233-1241.
    }

    \item Splines, interactions with splines, algebraic form of
    fitted model from \code{latex.rms}
    \enumerate{
      \item Knaus WA, Harrell FE, Lynn J, et al. (1995): The
      SUPPORT prognostic model: Objective estimates of
      survival for seriously ill hospitalized adults.  \emph{Annals
	of Internal Medicine} 122:191-203.
    }

    \item Splines, odds ratio chart from fitted model with
    nonlinear and interaction terms, use of \code{transcan} for
    imputation
    \enumerate{
      \item Lee KL, Woodlief LH, Topol EJ, Weaver WD, Betriu A.
      Col J, Simoons M, Aylward P, Van de Werf F, Califf RM.
      Predictors of 30-day mortality in the era of
      reperfusion for acute myocardial infarction: results
      from an international trial of 41,021 patients.
      \emph{Circulation} 1995;91:1659-1668.
    }

    \item Splines, external validation of logistic models,
    prediction rules using point tables
    \enumerate{
      \item Steyerberg EW, Hargrove YV, \emph{et al} (2001): Residual mass
      histology in testicular cancer: development and
      validation of a clinical prediction rule.  \emph{Stat in Med}
      2001;20:3847-3859.
      \item van Gorp MJ, Steyerberg EW, \emph{et al} (2003): Clinical
      prediction rule for 30-day mortality in Bjork-Shiley convexo-concave
      valve replacement.  \emph{J Clinical Epidemiology} 2003;56:1006-1012.
    }

    \item Model fitting, bootstrap validation, missing value
    imputation
    \enumerate{
      \item Krijnen P, van Jaarsveld BC, Steyerberg EW, Man in 't
      Veld AJ, Schalekamp, MADH, Habbema JDF (1998): A
      clinical prediction  rule for renal artery stenosis.
      \emph{Annals of Internal Medicine} 129:705-711.
    }

    \item Model fitting, splines, bootstrap validation, nomograms
    \enumerate{
      \item Kattan MW, Eastham JA, Stapleton AMF, Wheeler TM,
      Scardino PT.  A preoperative nomogram for disease
      recurrence following radical prostatectomy for
      prostate cancer.  \emph{J Natl Ca Inst} 1998;
      90(10):766-771.
      
      \item Kattan, MW, Wheeler TM, Scardino PT.  A
      postoperative nomogram for disease recurrence
      following radical prostatectomy for prostate
      cancer. \emph{J Clin Oncol} 1999; 17(5):1499-1507

      \item Kattan MW, Zelefsky MJ, Kupelian PA, Scardino PT, 
      Fuks Z, Leibel SA.  A pretreatment nomogram for
      predicting the outcome of three-dimensional
      conformal radiotherapy in prostate cancer.  
      \emph{J Clin Oncol} 2000; 18(19):3252-3259.
      
      \item Eastham JA, May R, Robertson JL, Sartor O, Kattan
      MW.  Development of a nomogram which predicts the
      probability of a positive prostate biopsy in men
      with an abnormal digital rectal examination and a
      prostate specific antigen between 0 and 4
      ng/ml. \emph{Urology}. (In press).
      
      \item Kattan MW, Heller G, Brennan MF.  A competing-risk
      nomogram fir sarcoma-specific death following local recurrence.
      \emph{Stat in Med} 2003; 22; 3515-3525.
    }

    \item Penalized maximum likelihood estimation, regression splines, web
    site to get predicted values
    \enumerate{
      \item Smits M, Dippel DWJ, Steyerberg EW, et al.  Predicting intracranial
      traumatic findings on computed tomography in patients with minor head
      injury: The CHIP prediction rule.  \emph{Ann Int Med} 2007; 146:397-405.
    }

    \item Nomogram with 2- and 5-year survival probability and median survival
    time (but watch out for the use of univariable screening)
    \enumerate{
      \item Clark TG, Stewart ME, Altman DG, Smyth JF.  A prognostic
      model for ovarian cancer.  \emph{Br J Cancer} 2001; 85:944-52.
    }

    \item Comprehensive example of parametric survival modeling
    with an extensive nomogram, time ratio chart, anova
    chart, survival curves generated using survplot,
    bootstrap calibration curve
    \enumerate{
      \item Teno JM, Harrell FE, Knaus WA, et al.  Prediction of
      survival for older hospitalized patients: The HELP
      survival model.  \emph{J Am Geriatrics Soc} 2000;
      48: S16-S24.
    }

    \item Model fitting, imputation, and several nomograms
    expressed in tabular form
    \enumerate{
      \item Hasdai D, Holmes DR, et al.  Cardiogenic shock complicating
      acute myocardial infarction: Predictors of death.
      \emph{Am Heart J} 1999; 138:21-31.
    }

    \item Ordinal logistic model with bootstrap calibration plot
    \enumerate{
      \item Wu AW, Yasui U, Alzola CF \emph{et al}.  Predicting functional
      status outcomes in hospitalized patients aged 80 years and
      older.  \emph{J Am Geriatric Society} 2000; 48:S6-S15.
    }

    \item Propensity modeling in evaluating medical diagnosis, anova
    dot chart
    \enumerate{
      \item Weiss JP, Gruver C, et al.  Ordering an echocardiogram 
      for evaluation of left ventricular function: Level
      of expertise necessary for efficient use. \emph{J Am Soc 
        Echocardiography} 2000; 13:124-130.
    }

    \item Simulations using rms to study the properties
    of various modeling strategies
    \enumerate{
      \item Steyerberg EW, Eijkemans MJC, Habbema JDF.  Stepwise selection
      in small data sets: A simulation study of bias in logistic
      regression analysis.  \emph{J Clin Epi} 1999; 52:935-942.

      \item Steyerberg WE, Eijekans MJC, Harrell FE, Habbema JDF.
      Prognostic modeling with logistic regression analysis: In
      search of a sensible strategy in small data sets.  \emph{Med
        Decision Making} 2001; 21:45-56.
    }

    \item Statistical methods and
    references related to rms, along with case studies
    which includes the rms code which produced the
    analyses
    \enumerate{
      \item Harrell FE, Lee KL, Mark DB (1996): Multivariable
      prognostic models: Issues in developing models,
      evaluating assumptions and adequacy, and measuring and
      reducing errors.  \emph{Stat in Med} 15:361-387.

      \item Harrell FE, Margolis PA, Gove S, Mason KE, Mulholland
      EK et al. (1998): Development of a clinical prediction
      model for an ordinal outcome: The World Health
      Organization ARI Multicentre Study of clinical signs
      and etiologic agents of pneumonia, sepsis, and
      meningitis in young infants. \emph{Stat in Med} 17:909-944.

      \item Bender R, Benner, A (2000): Calculating ordinal regression
      models in SAS and S-Plus.  \emph{Biometrical J} 42:677-699.
    }
  }
}

\section{Bug Reports}{
  The author is willing to help with problems.  Send
  E-mail to \email{f.harrell@vanderbilt.edu}.  To report bugs,
  please do the following:

  \enumerate{
    \item If the bug occurs when running a function on a fit
    object (e.g., \code{anova}), attach a \code{dump}'d text
    version of the fit object to your note.  If you
    used \code{datadist} but not until after the fit was
    created, also send the object created by
    \code{datadist}.  Example: \code{save(myfit,"/tmp/myfit.rda")} will create
    an R binary save file that can be attached to the E-mail.  
    \item If the bug occurs during a model fit (e.g., with
    \code{lrm, ols, psm, cph}), send the statement causing
    the error with a \code{save}'d version of the data
    frame used in the fit.  If this data frame is very
    large, reduce it to a small subset which still
    causes the error.
  }
}

\section{Copyright Notice}{
  GENERAL DISCLAIMER  This program is free software;
  you can redistribute it and/or modify it under the
  terms of the GNU General Public License as
  published by the Free Software Foundation; either
  version 2, or (at your option) any later version.

  This program is
  distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied
  warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.   In short: you may
  use this code any way you like, as long as you don't charge
  money for it, remove this notice, or hold anyone
  liable for its results.  Also, please acknowledge
  the source and communicate changes to the author.

  If this software is used is work presented for
  publication, kindly reference it using for example:
  Harrell FE (2009): rms: S functions for
  biostatistical/epidemiologic modeling, testing,
  estimation, validation, graphics, and prediction.
  Programs available from \url{biostat.mc.vanderbilt.edu/rms}.
  Be sure to reference other packages used as well as 
  \R itself.
}
 
\author{
  Frank E Harrell Jr\cr
  Professor of Biostatistics\cr
  Chair, Department of Biostatistics\cr
  Vanderbilt University School of Medicine\cr
  Nashville, Tennessee\cr
  \email{f.harrell@vanderbilt.edu}
}
\keyword{models}
\concept{overview}
