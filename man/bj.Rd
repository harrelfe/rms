\name{bj}
\alias{bj}
\alias{bj.fit}
\alias{residuals.bj}
\alias{print.bj}
\alias{validate.bj}
\alias{bjplot}
\title{
  Buckley-James Multiple Regression Model
}
\description{
  \code{bj} fits the Buckley-James distribution-free least squares multiple
  regression model to a possibly right-censored response variable.  
  This model reduces to ordinary least squares if
  there is no censoring.  By default, model fitting is done after
  taking logs of the response variable.
  \code{bj} uses the \code{rms} class
  for automatic \code{anova}, \code{fastbw}, \code{validate}, \code{Function}, \code{nomogram},
  \code{summary}, \code{plot}, \code{bootcov}, and other functions.  The \code{bootcov}
  function may be worth using with \code{bj} fits, as the properties of the
  Buckley-James covariance matrix estimator are not fully known for
  strange censoring patterns.

  The \code{residuals.bj} function exists mainly to compute 
  residuals and to censor them (i.e., return them as
  \code{Surv} objects) just as the original
  failure time variable was censored.  These residuals are useful for
  checking to see if the model also satisfies certain distributional assumptions.
  To get these residuals, the fit must have specified \code{y=TRUE}.

  The \code{bjplot} function is a special plotting function for objects
  created by \code{bj} with \code{x=TRUE, y=TRUE} in effect.  It produces three
  scatterplots for every covariate in the model: the first plots the
  original situation, where censored data are distingushed from
  non-censored data by a different plotting symbol. In the second plot,
  called a renovated plot, vertical lines show how censored data were
  changed by the procedure, and the third is equal to the second, but
  without vertical lines.  Imputed data are again distinguished from the
  non-censored by a different symbol.

  The \code{validate} method for \code{bj} validates the Somers' \code{Dxy} rank
  correlation between predicted and observed responses, accounting for censoring.

  The primary fitting function for \code{bj} is \code{bj.fit}, which does not
  allow missing data and expects a full design matrix as input.
}
\usage{
bj(formula=formula(data), data, subset, na.action=na.delete,
   link="log", control, method='fit', x=FALSE, y=FALSE, 
   time.inc)

\method{print}{bj}(x, digits=4, long=FALSE, coefs=TRUE, latex=FALSE,
title="Buckley-James Censored Data Regression", \dots)

\method{residuals}{bj}(object, type=c("censored","censored.normalized"),\dots)

bjplot(fit, which=1:dim(X)[[2]])

\method{validate}{bj}(fit, method="boot", B=40,
         bw=FALSE,rule="aic",type="residual",sls=.05,aics=0,
         force=NULL, estimates=TRUE, pr=FALSE,
		 tol=1e-7, rel.tolerance=1e-3, maxiter=15, \dots)

bj.fit(x, y, control)
}
\arguments{
  \item{formula}{
    an S statistical model formula. Interactions up to third order are
    supported. The left hand side must be a \code{Surv} object.
  }
  \item{data,subset,na.action}{the usual statistical model fitting arguments}
  \item{fit}{
    a fit created by \code{bj}, required for all functions except \code{bj}.
  }
  \item{x}{
    a design matrix with or without a first column of ones, to pass
    to \code{bj.fit}.  All models will have an intercept.  For
    \code{print.bj} is a result of \code{bj}.  For \code{bj}, set
    \code{x=TRUE} to include the design matrix in the fit object. 
  }
  \item{y}{
    a \code{Surv} object to pass to \code{bj.fit} as the two-column response 
    variable.  Only right censoring is allowed, and there need not be any
    censoring.  For \code{bj}, set \code{y} to \code{TRUE} to include the
    two-column response matrix, with the 
    event/censoring indicator in the second column.  The first column will
    be transformed according to \code{link}, and depending on
    \code{na.action}, rows with missing data in the predictors or the
    response will be deleted.
  }
  \item{link}{
    set to, for example, \code{"log"} (the default) to model the log of the
    response, or \code{"identity"} to model the untransformed response.
  }
  \item{control}{
    a list containing any or all of the following components: \code{iter.max}
    (maximum number of iterations allowed, default is 20),
    \code{eps} (convergence criterion: concergence is assumed when the ratio of
    sum of squared errors from one iteration to the next is between
    1-\code{eps} and 1+\code{eps}), \code{trace} (set to \code{TRUE} to monitor iterations), 
    \code{tol} (matrix singularity criterion, default is 1e-7), and 'max.cycle' 
    (in case of nonconvergence the program looks for a cycle that repeats itself, 
    default is 30).  
  }
  \item{method}{
    set to \code{"model.frame"} or \code{"model.matrix"} to return one of those
    objects rather than the model fit.
  }
  \item{time.inc}{
    setting for default time spacing.
    Default is 30 if time variable has \code{units="Day"}, 1 otherwise, unless
    maximum follow-up time \eqn{< 1}. Then max time/10 is used as \code{time.inc}.
    If \code{time.inc} is not given and max time/default \code{time.inc} is
    \eqn{> 25}, \code{time.inc} is increased.
  }
  \item{digits}{
    number of significant digits to print if not 4.
  }
  \item{long}{
    set to \code{TRUE} to print the correlation matrix for parameter estimates
  }
  \item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
	of model coefficients, standard errors, etc.  Specify \code{coefs=n}
	to print only the first \code{n} regression coefficients in the
	model.}
  \item{latex}{a logical value indicating whether information should be
	formatted as plain text or as LaTeX markup}
  \item{title}{a character string title to be passed to \code{prModFit}}
  \item{object}{the result of \code{bj}}
  \item{type}{
    type of residual desired.  Default is censored unnormalized residuals,
    defined as link(Y) - linear.predictors, where the
    link function was usually the log function.  You can specify
    \code{type="censored.normalized"} to divide the residuals by the estimate
    of \code{sigma}.
  }
  \item{which}{
    vector of integers or character strings naming elements of the design
    matrix (the names of the original predictors if they entered the model
    linearly) for which to have \code{bjplot} make plots of only the variables listed in \code{which} (names or numbers).
  }
  \item{B,bw,rule,sls,aics,force,estimates,pr,tol,rel.tolerance,maxiter}{see
	\code{\link{predab.resample}}} 
  \item{\dots}{
    ignored for \code{print}; passed through to
    \code{predab.resample} for \code{validate}
  }
}
\value{
  \code{bj} returns a fit object with similar information to what \code{survreg},
  \code{psm}, \code{cph} would store as 
  well as what \code{rms} stores and \code{units} and \code{time.inc}.
  \code{residuals.bj} returns a \code{Surv} object.  One of the components of the
  \code{fit} object produced by \code{bj} (and \code{bj.fit}) is a vector called
  \code{stats} which contains the following names elements: 
  \code{"Obs", "Events", "d.f.","error d.f.","sigma","g"}.  Here
  \code{sigma} is the estimate of the residual standard deviation.
  \code{g} is the \eqn{g}-index.  If the link function is \code{"log"},
  the \eqn{g}-index on the anti-log scale is also returned as \code{gr}.
}
\details{
  The program implements the algorithm as described in the original
  article by Buckley & James. Also, we have used the original Buckley &
  James prescription for computing variance/covariance estimator.  This
  is based on non-censored observations only and does not have any
  theoretical justification, but has been shown in simulation studies to
  behave well. Our experience confirms this view.  Convergence is rather
  slow with this method, so you may want to increase the number of
  iterations.  Our experience shows that often, in particular with high
  censoring, 100 iterations is not too many. Sometimes the method will
  not converge, but will instead enter a loop of repeating values (this
  is due to the discrete nature  
  of Kaplan and Meier estimator and usually happens with small sample sizes).
  The program will look for such a loop and return the average betas. It will also 
  issue a warning message and give the size of the cycle (usually less than 6).
}
\author{
  Janez Stare\cr
  Department of Biomedical Informatics\cr
  Ljubljana University\cr
  Ljubljana, Slovenia\cr
  \email{janez.stare@mf.uni-lj.si}


  Harald Heinzl\cr
  Department of Medical Computer Sciences\cr
  Vienna University\cr
  Vienna, Austria\cr
  \email{harald.heinzl@akh-wien.ac.at}


  Frank Harrell\cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
  \email{f.harrell@vanderbilt.edu}
}
\references{
  Buckley JJ, James IR. Linear regression with censored data. Biometrika 1979; 
  66:429--36.


  Miller RG, Halpern J. Regression with censored data. Biometrika 1982; 69: 
  521--31.


  James IR, Smith PJ. Consistency results for linear regression with censored 
  data. Ann Statist 1984; 12: 590--600.


  Lai TL, Ying Z. Large sample theory of a modified Buckley-James estimator for 
  regression analysis 
  with censored data. Ann Statist 1991; 19: 1370--402.


  Hillis SL. Residual plots for the censored data linear regression model.  Stat in Med 1995; 14: 2023--2036.


  Jin Z, Lin DY, Ying Z. On least-squares regression with censored data.  Biometrika 2006; 93:147--161.
}
\seealso{
  \code{\link{rms}}, \code{\link{psm}}, \code{\link{survreg}},
	\code{\link{cph}}, \code{\link[survival]{Surv}}, 
  \code{\link[Hmisc]{na.delete}},
  \code{\link[Hmisc]{na.detail.response}}, \code{\link{datadist}},
  \code{\link[Hmisc]{rcorr.cens}}, \code{\link{GiniMd}},
  \code{\link{prModFit}}, \code{\link{dxy.cens}}
}
\examples{
set.seed(1)
ftime  <- 10*rexp(200)
stroke <- ifelse(ftime > 10, 0, 1)
ftime  <- pmin(ftime, 10)
units(ftime) <- "Month"
age <- rnorm(200, 70, 10)
hospital <- factor(sample(c('a','b'),200,TRUE))
dd <- datadist(age, hospital)
options(datadist="dd")

f <- bj(Surv(ftime, stroke) ~ rcs(age,5) + hospital, x=TRUE, y=TRUE)
# add link="identity" to use a censored normal regression model instead
# of a lognormal one
anova(f)
fastbw(f)
validate(f, B=15)
plot(Predict(f, age, hospital))
# needs datadist since no explicit age,hosp.
coef(f)               # look at regression coefficients
coef(psm(Surv(ftime, stroke) ~ rcs(age,5) + hospital, dist='lognormal'))
                      # compare with coefficients from likelihood-based
                      # log-normal regression model
                      # use dist='gau' not under R 


r <- resid(f, 'censored.normalized')
survplot(npsurv(r ~ 1), conf='none') 
                      # plot Kaplan-Meier estimate of 
                      # survival function of standardized residuals
survplot(npsurv(r ~ cut2(age, g=2)), conf='none')  
                      # may desire both strata to be n(0,1)
options(datadist=NULL)
}
\keyword{models}
\keyword{survival}






