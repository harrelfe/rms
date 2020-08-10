\name{cph}
\alias{cph}
\alias{Survival.cph}
\alias{Quantile.cph}
\alias{Mean.cph}
\title{Cox Proportional Hazards Model and Extensions}
\description{
  Modification of Therneau's \code{coxph} function to fit the Cox model and
  its extension, the Andersen-Gill model. The latter allows for interval
  time-dependent covariables, time-dependent strata, and repeated events.
  The \code{Survival} method for an object created by \code{cph} returns an S
  function for computing estimates of the survival function.
  The \code{Quantile} method for \code{cph} returns an S function for computing
  quantiles of survival time (median, by default).
  The \code{Mean} method returns a function for computing the mean survival
  time.  This function issues a warning if the last follow-up time is uncensored,
  unless a restricted mean is explicitly requested.
}
\usage{
cph(formula = formula(data), data=environment(formula),
    weights, subset, na.action=na.delete, 
    method=c("efron","breslow","exact","model.frame","model.matrix"), 
    singular.ok=FALSE, robust=FALSE,
    model=FALSE, x=FALSE, y=FALSE, se.fit=FALSE,
    linear.predictors=TRUE, residuals=TRUE, nonames=FALSE,
    eps=1e-4, init, iter.max=10, tol=1e-9, surv=FALSE, time.inc,
    type=NULL, vartype=NULL, debug=FALSE, \dots)

\method{Survival}{cph}(object, \dots)
# Evaluate result as g(times, lp, stratum=1, type=c("step","polygon"))

\method{Quantile}{cph}(object, \dots)
# Evaluate like h(q, lp, stratum=1, type=c("step","polygon"))

\method{Mean}{cph}(object, method=c("exact","approximate"), type=c("step","polygon"),
          n=75, tmax, \dots)
# E.g. m(lp, stratum=1, type=c("step","polygon"), tmax, \dots)
}
\arguments{
  \item{formula}{
    an S formula object with a \code{Surv} object on the
          left-hand side.  The \code{terms} can specify any S model
          formula with up to third-order interactions.  The \code{strat} 
    function may appear in the terms, as a main effect or an interacting
    factor.  To stratify on both race and sex, you would include both
    terms \code{strat(race)} and \code{strat(sex)}.  Stratification
    factors may interact with non-stratification factors;
    not all stratification terms need interact with the same modeled
    factors.
  }
  \item{object}{
    an object created by \code{cph} with \code{surv=TRUE}
  }
  \item{data}{
    name of an S data frame containing all needed variables.  Omit this to use a
    data frame already in the S ``search list''.
  }
  \item{weights}{
    case weights
  }
  \item{subset}{
    an expression defining a subset of the observations to use in the fit.  The default
    is to use all observations.  Specify for example \code{age>50 & sex="male"} or
    \code{c(1:100,200:300)}
    respectively to use the observations satisfying a logical expression or those having
    row numbers in the given vector.
  }
  \item{na.action}{
    specifies an S function to handle missing data.  The default is the function \code{na.delete},
    which causes observations with any variable missing to be deleted.  The main difference
    between \code{na.delete} and the S-supplied function \code{na.omit} is that 
    \code{na.delete} makes a list
    of the number of observations that are missing on each variable in the model.
    The \code{na.action} is usally specified by e.g. \code{options(na.action="na.delete")}.
  }
  \item{method}{
    for \code{cph}, specifies a particular fitting method, \code{"model.frame"} instead to return the model frame
    of the predictor and response variables satisfying any subset or missing value
    checks, or \code{"model.matrix"} to return the expanded design matrix.
    The default is \code{"efron"}, to use Efron's likelihood for fitting the
    model.

    For \code{Mean.cph}, \code{method} is \code{"exact"} to use numerical
    integration of the 
    survival function at any linear predictor value to obtain a mean survival
    time.  Specify \code{method="approximate"} to use an approximate method that is
    slower when \code{Mean.cph} is executing but then is essentially instant
    thereafter.  For the approximate method, the area is computed for \code{n}
    points equally spaced between the min and max observed linear predictor
    values.  This calculation is done separately for each stratum.  Then the
    \code{n} pairs (X beta, area) are saved in the generated S function, and when
    this function is evaluated, the \code{approx} function is used to evaluate
    the mean for any given linear predictor values, using linear interpolation
    over the \code{n} X beta values.
  }
  \item{singular.ok}{
    If \code{TRUE}, the program will automatically skip over columns of the X matrix
    that are linear combinations of earlier columns.  In this case the
    coefficients for such columns will be NA, and the variance matrix will contain
    zeros.  For ancillary calculations, such as the linear predictor, the missing
    coefficients are treated as zeros.  The singularities will prevent many of
    the features of the \code{rms} library from working.
  }
  \item{robust}{
    if \code{TRUE} a robust variance estimate is returned.  Default is \code{TRUE} if the
    model includes a \code{cluster()} operative, \code{FALSE} otherwise.
  }
  \item{model}{
    default is \code{FALSE}(false).  Set to \code{TRUE} to return the model frame as element 
    \code{model} of the fit object.
  }
  \item{x}{
    default is \code{FALSE}.  Set to \code{TRUE} to return the expanded design matrix as element \code{x}
    (without intercept indicators) of the
    returned fit object.
  }
  \item{y}{
    default is \code{FALSE}.  Set to \code{TRUE} to return the vector of
    response values (\code{Surv} object) as element \code{y} of the
    fit.
  }
  \item{se.fit}{
    default is \code{FALSE}.  Set to \code{TRUE} to compute the estimated standard errors of
    the estimate of X beta and store them in element \code{se.fit}
    of the fit.  The predictors are first centered to their means
    before computing the standard errors.
  }
	\item{linear.predictors}{set to \code{FALSE} to omit
    \code{linear.predictors} vector from fit}
	\item{residuals}{set to \code{FALSE} to omit \code{residuals} vector
    from fit}
	\item{nonames}{set to \code{TRUE} to not set \code{names} attribute
    for \code{linear.predictors}, \code{residuals}, \code{se.fit}, and
    rows of design matrix}
  \item{eps}{
    convergence criterion - change in log likelihood.
  }
  \item{init}{
    vector of initial parameter estimates.  Defaults to all zeros.
    Special residuals can be obtained by setting some elements of \code{init}
    to MLEs and others to zero and specifying \code{iter.max=1}.
  }
  \item{iter.max}{
    maximum number of iterations to allow.  Set to \code{0} to obtain certain
    null-model residuals.
  }
  \item{tol}{
    tolerance for declaring singularity for matrix inversion (available
    only when survival5 or later package is in effect)
  }
  \item{surv}{
    set to \code{TRUE} to compute underlying survival estimates for each
    stratum, and to store these along with standard errors of log Lambda(t),
    \code{maxtime} (maximum observed survival or censoring time),
    and \code{surv.summary} in the returned object.  Set \code{surv="summary"}
    to only compute and store \code{surv.summary}, not survival estimates
    at each unique uncensored failure time. If you specify \code{x=TRUE}
    and \code{y=TRUE}, 
    you can obtain predicted survival later, with accurate confidence
    intervals for any set of predictor values. The standard error information
    stored as a result of \code{surv=TRUE} are only accurate at the mean of all
    predictors. If the model has no covariables, these are of course OK.
    The main reason for using \code{surv} is to greatly speed up the computation
    of predicted survival probabilities as a function of the covariables,
    when accurate confidence intervals are not needed.
  }
  \item{time.inc}{
    time increment used in deriving \code{surv.summary}.  Survival,
    number at risk, and standard error will be stored for 
    \code{t=0, time.inc, 2 time.inc, \dots, maxtime},
    where \code{maxtime} is the maximum survival time over all strata.
    \code{time.inc} is also used in constructing the time axis in the
    \code{survplot} function (see below).  The default value for
    \code{time.inc} is 30 if \code{units(ftime) = "Day"} or no \code{units}
    attribute has been attached to the survival time variable.  If
    \code{units(ftime)} is a word other than \code{"Day"}, the default
    for \code{time.inc} is 1 when it is omitted, unless \code{maxtime<1}, then
    \code{maxtime/10} is used as \code{time.inc}.  If \code{time.inc} is not given and
    \code{maxtime/ default time.inc} > 25, \code{time.inc} is increased.
  }
  \item{type}{
    (for \code{cph}) applies if \code{surv} is \code{TRUE} or \code{"summary"}. 
    If \code{type} is omitted, the method consistent with \code{method} is used.
    See \code{survfit.coxph} (under \code{survfit}) or \code{survfit.cph} for details and for the
    definitions of values of \code{type}

    For \code{Survival, Quantile, Mean} set to \code{"polygon"} to use linear 
    interpolation instead of the usual step function.  For \code{Mean}, the default
    of \code{step} will yield the sample mean in the case of no censoring and no
    covariables, if \code{type="kaplan-meier"} was specified to \code{cph}.
    For \code{method="exact"}, the value of \code{type} is passed to the
    generated function, and it can be overridden when that function is
    actually invoked. For \code{method="approximate"}, \code{Mean.cph}
    generates the function different ways according to \code{type}, and this
    cannot be changed when the function is actually invoked.
  }
  \item{vartype}{see \code{survfit.coxph}}
	\item{debug}{set to \code{TRUE} to print debugging information related
    to model matrix construction.  You can also use \code{options(debug=TRUE)}.}
  \item{\dots}{
    other arguments passed to \code{coxph.fit} from \code{cph}.  Ignored by
    other functions.
  }
  \item{times}{
    a scalar or vector of times at which to evaluate the survival estimates
  }
  \item{lp}{
    a scalar or vector of linear predictors (including the centering constant)
    at which to evaluate the survival estimates
  }
  \item{stratum}{
    a scalar stratum number or name (e.g., \code{"sex=male"}) to use in getting
    survival probabilities
  }
  \item{q}{
    a scalar quantile or a vector of quantiles to compute
  }
  \item{n}{
    the number of points at which to evaluate the mean survival time, for
    \code{method="approximate"} in \code{Mean.cph}.
  }
  \item{tmax}{
    For \code{Mean.cph}, the default is to compute the overall mean (and produce
    a warning message if there is censoring at the end of follow-up).
    To compute a restricted mean life length, specify the truncation point as \code{tmax}.
    For \code{method="exact"}, \code{tmax} is passed to the generated function and it
    may be overridden when that function is invoked.  For \code{method="approximate"},
    \code{tmax} must be specified at the time that \code{Mean.cph} is run.
}}
\value{
  For \code{Survival}, \code{Quantile}, or \code{Mean}, an S function is returned.  Otherwise,
  in addition to what is listed below, formula/design information and
  the components 
  \code{maxtime, time.inc, units, model, x, y, se.fit} are stored, the last 5 
  depending on the settings of options by the same names.
  The vectors or matrix stored if \code{y=TRUE} or \code{x=TRUE} have rows deleted according to \code{subset} and
  to missing data, and have names or row names that come from the
  data frame used as input data.

  \item{n}{
    table with one row per stratum containing number of censored and uncensored observations
  }
  \item{coef}{
    vector of regression coefficients
  }
  \item{stats}{
    vector containing the named elements \code{Obs}, \code{Events}, \code{Model L.R.}, \code{d.f.},
    \code{P}, \code{Score}, \code{Score P}, \code{R2}, Somers'
    \code{Dxy}, \code{g}-index, 
    and \code{gr}, the \code{g}-index on the hazard ratio scale.
    \code{R2} is the Nagelkerke R-squared, with division by the maximum
    attainable R-squared.
  }
  \item{var}{
    variance/covariance matrix of coefficients
  }
  \item{linear.predictors}{
    values of predicted X beta for observations used in fit, normalized
    to have overall mean zero, then having any offsets added
  }
  \item{resid}{
    martingale residuals
  }
  \item{loglik}{
    log likelihood at initial and final parameter values
  }
  \item{score}{
    value of score statistic at initial values of parameters
  }
  \item{times}{
    lists of times (if \code{surv="T"})
  }
  \item{surv}{
    lists of underlying survival probability estimates
  }
  \item{std.err}{
    lists of standard errors of estimate log-log survival
  }
  \item{surv.summary}{
    a 3 dimensional array if \code{surv=TRUE}.  
    The first dimension is time ranging from 0 to
    \code{maxtime} by \code{time.inc}.  The second dimension refers to strata.
    The third dimension contains the time-oriented matrix with
    \code{Survival, n.risk} (number of subjects at risk), 
    and \code{std.err} (standard error of log-log
    survival). 
  }
  \item{center}{
    centering constant, equal to overall mean of X beta.
}}
\details{
  If there is any strata by covariable interaction in the model such that
  the mean X beta varies greatly over strata, \code{method="approximate"} may
  not yield very accurate estimates of the mean in \code{Mean.cph}.


  For \code{method="approximate"} if you ask for an estimate of the mean for
  a linear predictor value that was outside the range of linear predictors
  stored with the fit, the mean for that observation will be \code{NA}.
}
\author{
  Frank Harrell\cr
  Department of Biostatistics, Vanderbilt University\cr
  \email{fh@fharrell.com}
}
\seealso{
  \code{\link[survival]{coxph}}, \code{\link[survival]{survival-internal}},
  \code{\link[survival]{Surv}}, \code{\link{residuals.cph}},
  \code{\link[survival]{cox.zph}}, \code{\link{survfit.cph}},
  \code{\link{survest.cph}}, \code{\link[survival]{survfit.coxph}}, 
  \code{\link{survplot}}, \code{\link{datadist}},
  \code{\link{rms}}, \code{\link{rms.trans}}, \code{\link{anova.rms}},
  \code{\link{summary.rms}}, \code{\link{Predict}}, 
  \code{\link{fastbw}}, \code{\link{validate}}, \code{\link{calibrate}},
  \code{\link{plot.Predict}}, \code{\link{ggplot.Predict}},
	\code{\link{specs.rms}}, \code{\link{lrm}}, \code{\link{which.influence}},
  \code{\link[Hmisc]{na.delete}},
  \code{\link[Hmisc]{na.detail.response}},  \code{\link{print.cph}},
  \code{\link{latex.cph}}, \code{\link{vif}}, \code{\link{ie.setup}},
  \code{\link{GiniMd}}, \code{\link{dxy.cens}},
  \code{\link[survival:survConcordance]{survConcordance}}
}
\examples{
# Simulate data from a population model in which the log hazard
# function is linear in age and there is no age x sex interaction

n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
              rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
dd <- datadist(age, sex)
options(datadist='dd')
S <- Surv(dt,e)

f <- cph(S ~ rcs(age,4) + sex, x=TRUE, y=TRUE)
cox.zph(f, "rank")             # tests of PH
anova(f)
ggplot(Predict(f, age, sex)) # plot age effect, 2 curves for 2 sexes
survplot(f, sex)             # time on x-axis, curves for x2
res <- resid(f, "scaledsch")
time <- as.numeric(dimnames(res)[[1]])
z <- loess(res[,4] ~ time, span=0.50)   # residuals for sex
plot(time, fitted(z))
lines(supsmu(time, res[,4]),lty=2)
plot(cox.zph(f,"identity"))    #Easier approach for last few lines
# latex(f)


f <- cph(S ~ age + strat(sex), surv=TRUE)
g <- Survival(f)   # g is a function
g(seq(.1,1,by=.1), stratum="sex=Male", type="poly") #could use stratum=2
med <- Quantile(f)
plot(Predict(f, age, fun=function(x) med(lp=x)))  #plot median survival

# Fit a model that is quadratic in age, interacting with sex as strata
# Compare standard errors of linear predictor values with those from
# coxph
# Use more stringent convergence criteria to match with coxph

f <- cph(S ~ pol(age,2)*strat(sex), x=TRUE, eps=1e-9, iter.max=20)
coef(f)
se <- predict(f, se.fit=TRUE)$se.fit
require(lattice)
xyplot(se ~ age | sex, main='From cph')
a <- c(30,50,70)
comb <- data.frame(age=rep(a, each=2),
                   sex=rep(levels(sex), 3))

p <- predict(f, comb, se.fit=TRUE)
comb$yhat  <- p$linear.predictors
comb$se    <- p$se.fit
z <- qnorm(.975)
comb$lower <- p$linear.predictors - z*p$se.fit
comb$upper <- p$linear.predictors + z*p$se.fit
comb

age2 <- age^2
f2 <- coxph(S ~ (age + age2)*strata(sex))
coef(f2)
se <- predict(f2, se.fit=TRUE)$se.fit
xyplot(se ~ age | sex, main='From coxph')
comb <- data.frame(age=rep(a, each=2), age2=rep(a, each=2)^2,
                   sex=rep(levels(sex), 3))
p <- predict(f2, newdata=comb, se.fit=TRUE)
comb$yhat <- p$fit
comb$se   <- p$se.fit
comb$lower <- p$fit - z*p$se.fit
comb$upper <- p$fit + z*p$se.fit
comb


# g <- cph(Surv(hospital.charges) ~ age, surv=TRUE)
# Cox model very useful for analyzing highly skewed data, censored or not
# m <- Mean(g)
# m(0)                           # Predicted mean charge for reference age


#Fit a time-dependent covariable representing the instantaneous effect
#of an intervening non-fatal event
rm(age)
set.seed(121)
dframe <- data.frame(failure.time=1:10, event=rep(0:1,5),
                     ie.time=c(NA,1.5,2.5,NA,3,4,NA,5,5,5), 
                     age=sample(40:80,10,rep=TRUE))
z <- ie.setup(dframe$failure.time, dframe$event, dframe$ie.time)
S <- z$S
ie.status <- z$ie.status
attach(dframe[z$subs,])    # replicates all variables

f <- cph(S ~ age + ie.status, x=TRUE, y=TRUE)  
#Must use x=TRUE,y=TRUE to get survival curves with time-dep. covariables


#Get estimated survival curve for a 50-year old who has an intervening
#non-fatal event at 5 days
new <- data.frame(S=Surv(c(0,5), c(5,999), c(FALSE,FALSE)), age=rep(50,2),
                  ie.status=c(0,1))
g <- survfit(f, new)
plot(c(0,g$time), c(1,g$surv[,2]), type='s', 
     xlab='Days', ylab='Survival Prob.')
# Not certain about what columns represent in g$surv for survival5
# but appears to be for different ie.status
#or:
#g <- survest(f, new)
#plot(g$time, g$surv, type='s', xlab='Days', ylab='Survival Prob.')


#Compare with estimates when there is no intervening event
new2 <- data.frame(S=Surv(c(0,5), c(5, 999), c(FALSE,FALSE)), age=rep(50,2),
                   ie.status=c(0,0))
g2 <- survfit(f, new2)
lines(c(0,g2$time), c(1,g2$surv[,2]), type='s', lty=2)
#or:
#g2 <- survest(f, new2)
#lines(g2$time, g2$surv, type='s', lty=2)
detach("dframe[z$subs, ]")
options(datadist=NULL)
}
\keyword{survival}
\keyword{models}
\keyword{nonparametric}
