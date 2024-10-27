\name{psm}
\alias{psm}
\alias{print.psm}
\alias{Hazard}
\alias{Survival}
\alias{Hazard.psm}
\alias{Mean.psm}
\alias{Quantile.psm}
\alias{Survival.psm}
\alias{residuals.psm}
\alias{lines.residuals.psm.censored.normalized}
\alias{survplot.residuals.psm.censored.normalized}
\title{Parametric Survival Model}

\description{
  \code{psm} is a modification of Therneau's \code{survreg} function for
  fitting the accelerated failure time family of parametric survival
  models.  \code{psm} uses the \code{rms} class for automatic
  \code{anova}, \code{fastbw}, \code{calibrate}, \code{validate}, and
  other functions.  \code{Hazard.psm}, \code{Survival.psm},
  \code{Quantile.psm}, and \code{Mean.psm} create S functions that
  evaluate the hazard, survival, quantile, and mean (expected value)
  functions analytically, as functions of time or probabilities and the
  linear predictor values.  The Nagelkerke R^2 and and adjusted
  Maddala-Cox-Snell R^2 are computed.  For the latter the notation is
  R2(p,m) where p is the number of regression coefficients being
  adjusted for and m is the effective sample size (number of uncensored
  observations).  See \code{\link[Hmisc]{R2Measures}} for more information.

	For the \code{print} method, format of output is controlled by the
	user previously running \code{options(prType="lang")} where
	\code{lang} is \code{"plain"} (the default), \code{"latex"}, or
	\code{"html"}. 


  The \code{residuals.psm} function exists mainly to compute normalized
  (standardized) residuals and to censor them (i.e., return them as
  \code{Surv} objects) just as the original failure time variable was
  censored.  These residuals are useful for checking the underlying
  distributional assumption (see the examples).  To get these residuals,
  the fit must have specified \code{y=TRUE}.  A \code{lines} method for these
  residuals automatically draws a curve with the assumed standardized
  survival distribution.  A \code{survplot} method runs the standardized
  censored residuals through \code{npsurv} to get Kaplan-Meier estimates,
  with optional stratification (automatically grouping a continuous
  variable into quantiles) and then through \code{survplot.npsurv} to plot
  them.  Then \code{lines} is invoked to show the theoretical curve.  Other
  types of residuals are computed by \code{residuals} using
  \code{residuals.survreg}.
}
\usage{
psm(formula,
    data=environment(formula), weights,
    subset, na.action=na.delete, dist="weibull",
    init=NULL, scale=0, 
    control=survreg.control(),
    parms=NULL, 
    model=FALSE, x=FALSE, y=TRUE, time.inc, \dots)

\method{print}{psm}(x, correlation=FALSE, digits=4, r2=c(0,2,4), coefs=TRUE,
pg=FALSE, title, \dots) 

Hazard(object, \dots)
\method{Hazard}{psm}(object, \dots)   # for psm fit
# E.g. lambda <- Hazard(fit)

Survival(object, \dots)
\method{Survival}{psm}(object, \dots) # for psm
# E.g. survival <- Survival(fit)

\method{Quantile}{psm}(object, \dots) # for psm
# E.g. quantsurv <- Quantile(fit)

\method{Mean}{psm}(object, \dots)     # for psm
# E.g. meant   <- Mean(fit)

# lambda(times, lp)   # get hazard function at t=times, xbeta=lp
# survival(times, lp) # survival function at t=times, lp
# quantsurv(q, lp)    # quantiles of survival time
# meant(lp)           # mean survival time

\method{residuals}{psm}(object, type=c("censored.normalized",
"response", "deviance", "dfbeta",
"dfbetas", "working", "ldcase", "ldresp", "ldshape", "matrix", "score"), \dots)

\method{survplot}{residuals.psm.censored.normalized}(fit, x, g=4, col, main, \dots)

\method{lines}{residuals.psm.censored.normalized}(x, n=100, lty=1, xlim,
lwd=3, \dots)
# for type="censored.normalized"
}
\arguments{
  \item{formula}{
    an S statistical model formula. Interactions up to third order are
    supported. The left hand side must be a \code{Surv} object.
  }
  \item{object}{a fit created by \code{psm}.  For \code{survplot} with
    residuals from \code{psm}, \code{object} is the result of
    \code{residuals.psm}.
  }
  \item{fit}{a fit created by \code{psm}}
  \item{data,subset,weights,dist,scale,init,na.action,control}{see \code{survreg}.}
  \item{parms}{a list of fixed parameters.  For the \eqn{t}-distribution
    this is the degrees of freedom; most of the distributions have no
    parameters.}
  \item{model}{
    set to \code{TRUE} to include the model frame in the returned object
  }
  \item{x}{
    set to \code{TRUE} to include the design matrix in the object produced
    by \code{psm}.  For the \code{survplot} method, \code{x} is an optional
    stratification variable (character, numeric, or categorical).  For
    \code{lines.residuals.psm.censored.normalized}, \code{x} is the result
    of \code{residuals.psm}.  For \code{print} it is the result of \code{psm}.
  }
  \item{y}{
    set to \code{TRUE} to include the \code{Surv()} matrix
  }
  \item{time.inc}{
    setting for default time spacing. Used in constructing time axis
    in \code{survplot}, and also in make confidence bars. Default is 30
    if time variable has \code{units="Day"}, 1 otherwise, unless
    maximum follow-up time \eqn{< 1}. Then max time/10 is used as \code{time.inc}.
    If \code{time.inc} is not given and max time/default \code{time.inc} is
    \eqn{> 25}, \code{time.inc} is increased.
  }
  \item{correlation}{set to \code{TRUE} to print the correlation matrix
    for parameter estimates}
  \item{digits}{number of places to print to the right of the decimal
		point}
	\item{r2}{vector of integers specifying which R^2 measures to print,
		with 0 for Nagelkerke R^2 and 1:4 corresponding to the 4 measures
		computed by \code{\link[Hmisc]{R2Measures}}.  Default is to print
		Nagelkerke (labeled R2) and second and fourth \code{R2Measures}
		which are the measures adjusted for the number of predictors, first
		for the raw sample size then for the effective sample size, which
		here is the number of uncensored observations.}
  \item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
	of model coefficients, standard errors, etc.  Specify \code{coefs=n}
	to print only the first \code{n} regression coefficients in the
	model.}
  \item{pg}{set to \code{TRUE} to print g-indexes}
  \item{title}{a character string title to be passed to \code{prModFit}}
  \item{\dots}{
    other arguments to fitting routines, or to pass to \code{survplot} from
    \cr
    \code{survplot.residuals.psm.censored.normalized}.  Passed to the
	generic \code{lines} function for \code{lines}.}
  \item{times}{
    a scalar or vector of times for which to evaluate survival probability
    or hazard
  }
  \item{lp}{
    a scalar or vector of linear predictor values at which to evaluate
    survival probability or hazard.  If both \code{times} and \code{lp} are
    vectors, they must be of the same length.
  }
  \item{q}{
    a scalar or vector of probabilities.  The default is .5, so just the
    median survival time is returned.  If \code{q} and \code{lp} are both vectors,
    a matrix of quantiles is returned, with rows corresponding to \code{lp}
    and columns to \code{q}.
  }
  \item{type}{
    type of residual desired.  Default is censored normalized residuals,
    defined as (link(Y) - linear.predictors)/scale parameter, where the
    link function was usually the log function.  See \code{survreg} for other
    types.  \code{type="score"} returns the score residual matrix.
  }
  \item{n}{
    number of points to evaluate theoretical standardized survival
    function for 
    \cr
    \code{lines.residuals.psm.censored.normalized}
  }
  \item{lty}{
    line type for \code{lines}, default is 1
  }
  \item{xlim}{
    range of times (or transformed times) for which to evaluate the standardized
    survival function.  Default is range in normalized residuals.
  }
  \item{lwd}{
    line width for theoretical distribution, default is 3
  }
  \item{g}{
    number of quantile groups to use for stratifying continuous variables
    having more than 5 levels
  }
  \item{col}{
    vector of colors for \code{survplot} method, corresponding to levels of \code{x}
    (must be a scalar if there is no \code{x})
  }
  \item{main}{
    main plot title for \code{survplot}.  If omitted, is the name or label of
    \code{x} if \code{x} is given.  Use \code{main=""} to suppress a title when you
    specify \code{x}.
}}
\value{
  \code{psm} returns a fit object with all the information \code{survreg} would store as 
  well as what \code{rms} stores and \code{units} and \code{time.inc}.
  \code{Hazard}, \code{Survival}, and \code{Quantile} return S-functions.
  \code{residuals.psm} with \code{type="censored.normalized"} returns a
  \code{Surv} object which has a special attribute \code{"theoretical"}
  which is used by the \code{lines} 
  routine.  This is the assumed standardized survival function as a function
  of time or transformed time.
}
\details{
  The object \code{survreg.distributions} contains definitions of properties
  of the various survival distributions. 
  \cr
  \code{psm} does not trap singularity errors due to the way \code{survreg.fit}
  does matrix inversion.  It will trap non-convergence (thus returning
  \code{fit$fail=TRUE}) if you give the argument \code{failure=2} inside the
  \code{control} list which is passed to \code{survreg.fit}.  For example, use
  \code{f <- psm(S ~ x, control=list(failure=2, maxiter=20))} to allow up to
  20 iterations and to set \code{f$fail=TRUE} in case of non-convergence.
  This is especially useful in simulation work.
}
\author{
  Frank Harrell\cr
  Department of Biostatistics\cr
  Vanderbilt University
  \cr
  \email{fh@fharrell.com}
}
\seealso{
  \code{\link{rms}}, \code{\link[survival]{survreg}},
  \code{\link[survival]{residuals.survreg}}, \code{\link[survival]{survreg.object}},  
  \code{\link[survival]{survreg.distributions}},
  \code{\link{pphsm}}, \code{\link{survplot}}, \code{\link{survest}},
  \code{\link[survival]{Surv}}, 
  \code{\link[Hmisc]{na.delete}},
  \code{\link[Hmisc]{na.detail.response}}, \code{\link{datadist}},
  \code{\link{latex.psm}}, \code{\link[Hmisc]{GiniMd}}, \code{\link{prModFit}},
	\code{\link{ggplot.Predict}}, \code{\link{plot.Predict}},
	\code{\link[Hmisc]{R2Measures}}
}
\examples{
require(survival)
n <- 400
set.seed(1)
age <- rnorm(n, 50, 12)
sex <- factor(sample(c('Female','Male'),n,TRUE))
dd <- datadist(age,sex)
options(datadist='dd')
# Population hazard function:
h <- .02*exp(.06*(age-50)+.8*(sex=='Female'))
d.time <- -log(runif(n))/h
cens <- 15*runif(n)
death <- ifelse(d.time <= cens,1,0)
d.time <- pmin(d.time, cens)

f <- psm(Surv(d.time,death) ~ sex*pol(age,2), 
         dist='lognormal')
# Log-normal model is a bad fit for proportional hazards data
print(f, r2=0:4, pg=TRUE)

anova(f)
fastbw(f)  # if deletes sex while keeping age*sex ignore the result
f <- update(f, x=TRUE,y=TRUE)       # so can validate, compute certain resids
validate(f, B=10)      # ordinarily use B=300 or more
plot(Predict(f, age, sex))   # needs datadist since no explicit age, hosp.
# Could have used ggplot(Predict(...))
survplot(f, age=c(20,60))     # needs datadist since hospital not set here
# latex(f)


S <- Survival(f)
plot(f$linear.predictors, S(6, f$linear.predictors),
     xlab=expression(X*hat(beta)),
     ylab=expression(S(6,X*hat(beta))))
# plots 6-month survival as a function of linear predictor (X*Beta hat)


times <- seq(0,24,by=.25)
plot(times, S(times,0), type='l')   # plots survival curve at X*Beta hat=0
lam <- Hazard(f)
plot(times, lam(times,0), type='l') # similarly for hazard function


med <- Quantile(f)        # new function defaults to computing median only
lp <- seq(-3, 5, by=.1)
plot(lp, med(lp=lp), ylab="Median Survival Time")
med(c(.25,.5), f$linear.predictors)
                          # prints matrix with 2 columns


# fit a model with no predictors
f <- psm(Surv(d.time,death) ~ 1, dist="weibull")
f
pphsm(f)          # print proportional hazards form
g <- survest(f)
plot(g$time, g$surv, xlab='Time', type='l',
     ylab=expression(S(t)))


f <- psm(Surv(d.time,death) ~ age, 
         dist="loglogistic", y=TRUE)
r <- resid(f, 'cens') # note abbreviation
survplot(npsurv(r ~ 1), conf='none') 
                      # plot Kaplan-Meier estimate of 
                      # survival function of standardized residuals
survplot(npsurv(r ~ cut2(age, g=2)), conf='none')  
                      # both strata should be n(0,1)
lines(r)              # add theoretical survival function
#More simply:
survplot(r, age, g=2)

options(datadist=NULL)
}
\keyword{models}
\keyword{survival}
