\name{validate.cph}
\alias{validate.cph}
\alias{validate.psm}
\title{Validation of a Fitted Cox or Parametric Survival Model's Indexes
  of Fit}
\description{
  This is the version of the \code{validate} function specific to models
  fitted with \code{cph} or \code{psm}. 
}
\usage{
# fit <- cph(formula=Surv(ftime,event) ~ terms, x=TRUE, y=TRUE, \dots)
\method{validate}{cph}(fit, method="boot", B=40, bw=FALSE, rule="aic", type="residual",
        sls=.05, aics=0, force=NULL, pr=FALSE, dxy=FALSE, u, tol=1e-9, \dots)

\method{validate}{psm}(fit, method="boot",B=40,
        bw=FALSE, rule="aic", type="residual", sls=.05, aics=0,
        force=NULL,
        pr=FALSE,
        dxy=FALSE, tol=1e-12, rel.tolerance=1e-5, maxiter=15, \dots)
}
\arguments{
  \item{fit}{
    a fit derived \code{cph}. The options \code{x=TRUE} and \code{y=TRUE}
    must have been specified. If the model contains any stratification factors
    and dxy=TRUE,
    the options \code{surv=TRUE} and \code{time.inc=u} must also have been given,
    where \code{u} is the same value of \code{u} given to \code{validate}.
  }
  \item{method}{see \code{\link{validate}}}
  \item{B}{
    number of repetitions.  For \code{method="crossvalidation"}, is the
    number of groups of omitted observations.
  }
  \item{rel.tolerance,maxiter,bw}{
    \code{TRUE} to do fast step-down using the \code{fastbw} function,
    for both the overall model and for each repetition. \code{fastbw}
    keeps parameters together that represent the same factor.
  }
  \item{rule}{
    Applies if \code{bw=TRUE}.  \code{"aic"} to use Akaike's information criterion as a
    stopping rule (i.e., a factor is deleted if the \eqn{\chi^2}{chi-square} falls below
    twice its degrees of freedom), or \code{"p"} to use \eqn{P}-values.
  }
  \item{type}{
    \code{"residual"} or \code{"individual"} - stopping rule is for individual factors or
    for the residual \eqn{\chi^2}{chi-square} for all variables deleted
  }
  \item{sls}{
    significance level for a factor to be kept in a model, or for judging the
    residual \eqn{\chi^2}{chi-square}.
  }
  \item{aics}{
    cutoff on AIC when \code{rule="aic"}.
  }
  \item{force}{see \code{\link{fastbw}}}
  \item{pr}{
    \code{TRUE} to print results of each repetition
  }
  \item{tol,\dots}{see \code{\link{validate}} or \code{\link{predab.resample}}}
  \item{dxy}{
    set to \code{TRUE} to validate Somers' \eqn{D_{xy}}{Dxy}  using
    \code{rcorr.cens}, which takes longer.
  }
  \item{u}{
    must be specified if the model has any stratification factors and \code{dxy=TRUE}.
    In that case, strata are not included in \eqn{X\beta}{X beta} and the
    survival curves may cross.  Predictions at time \code{t=u} are
    correlated with observed survival times.  Does not apply to
    \code{validate.psm}.
  }
}
\details{
  Statistics validated include the Nagelkerke \eqn{R^2}, 
  \eqn{D_{xy}}{Dxy}, slope shrinkage,  the discrimination index \eqn{D}
  [(model L.R. \eqn{\chi^2}{chi-square} - 1)/L], the unreliability index
  \eqn{U} = (difference in -2 log likelihood between uncalibrated
  \eqn{X\beta}{X beta} and  
  \eqn{X\beta}{X beta} with overall slope calibrated to test sample) / L,
  and the overall quality index \eqn{Q = D - U}.  \eqn{g} is the
        \eqn{g}-index on the log relative hazard (linear predictor) scale.
  L is -2 log likelihood with beta=0.  The "corrected" slope
  can be thought of as shrinkage factor that takes into account overfitting.
  See \code{predab.resample} for the list of resampling methods.
}
\value{
  matrix with rows corresponding to \eqn{D_{xy}}{Dxy}, Slope, \eqn{D},
  \eqn{U}, and \eqn{Q}, and columns for the original index, resample estimates, 
  indexes applied to whole or omitted sample using model derived from
  resample, average optimism, corrected index, and number of successful
  resamples.\cr

  The values corresponting to the row \eqn{D_{xy}}{Dxy} are equal to \eqn{2 *
    (C - 0.5)} where C is the C-index or concordance probability.
  If the user is correlating the linear predictor (predicted log hazard)
  with survival time and she wishes to get the more usual correlation
  using predicted survival time or predicted survival probability,
  \eqn{D_{xy}}{Dxy} should be negated.
  
}
\section{Side Effects}{
  prints a summary, and optionally statistics for each re-fit (if \code{pr=TRUE})
}
\author{
  Frank Harrell\cr
  Department of Biostatistics, Vanderbilt University\cr
  f.harrell@vanderbilt.edu
}
\seealso{
  \code{\link{validate}}, \code{\link{predab.resample}},
  \code{\link{fastbw}}, \code{\link{rms}}, \code{\link{rms.trans}},
  \code{\link{calibrate}}, \code{\link[Hmisc]{rcorr.cens}},
  \code{\link{cph}}, \code{\link[survival]{survival-internal}},
  \code{\link{gIndex}}
}
\examples{
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, TRUE))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
S <- Surv(dt,e)

f <- cph(S ~ age*sex, x=TRUE, y=TRUE)
# Validate full model fit
validate(f, B=10)               # normally B=150

# Validate a model with stratification.  Dxy is the only
# discrimination measure for such models, by Dxy requires
# one to choose a single time at which to predict S(t|X)
f <- cph(S ~ rcs(age)*strat(sex), 
         x=TRUE, y=TRUE, surv=TRUE, time.inc=2)
validate(f, dxy=TRUE, u=2, B=10)   # normally B=150
# Note u=time.inc
}
\keyword{models}
\keyword{regression}
\keyword{survival}
\concept{model validation}
\concept{predictive accuracy}
\concept{bootstrap}
