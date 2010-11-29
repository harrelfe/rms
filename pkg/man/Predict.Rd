\name{Predict}
\alias{Predict}
\alias{print.Predict}
\alias{rbind.Predict}
\title{Compute Predicted Values and Confidence Limits}
\description{
\code{Predict} allows the user to easily specify which predictors are to
vary.  When the vector of values over which a predictor should vary is
not specified, the
range will be all levels of a categorical predictor or equally-spaced
points between the \code{\link{datadist}} \code{"Low:prediction"} and
\code{"High:prediction"} values for the variable (\code{datadist} by
default uses the 10th smallest and 10th largest predictor values in the
dataset).  Predicted values are 
the linear predictor (X beta), a user-specified transformation of that
scale, or estimated probability of surviving past a fixed single time
point given the linear predictor.  \code{Predict} is usually used for
plotting predicted values but there is also a \code{print} method.

There is a \code{plot} method for \code{Predict} objects that makes it
easy to show predicted values and confidence bands.

The \code{rbind} method for \code{Predict} objects allows you to create
separate sets of predictions under different situations and to combine
them into one set for feeding to \code{plot.Predict}.  For example you
might want to plot confidence intervals for means and for individuals
using \code{ols}, and have the two types of confidence bands be
superposed onto one plot or placed into two panels.  Another use for
\code{rbind} is to combine predictions from quantile regression models
that predicted three different quantiles.
}
\usage{
Predict(x, ..., fun,
        type = c("predictions", "model.frame", "x"),
        np = 200, conf.int = 0.95,
        conf.type = c("mean", "individual"),
        adj.zero = FALSE, ref.zero = FALSE,
        non.slopes, time = NULL, loglog = FALSE, digits=4, name, factors=NULL)

\method{print}{Predict}(x, \dots)

\method{rbind}{Predict}(\dots, rename)
}
\arguments{
  \item{x}{
	an \code{rms} fit object, or for \code{print} the result of \code{Predict}.
	\code{options(datadist="d")} must have been specified (where
	\code{d} was created by \code{datadist}), or 
	it must have been in effect when the the model was fitted.}
  \item{\dots}{
	One or more variables to vary, or single-valued adjustment values.
	Specify a variable name without an equal sign to use the default
	display range, or any range 
	you choose (e.g. \code{seq(0,100,by=2),c(2,3,7,14)}). 
	The default list of values for which predictions are made
	is taken as the list of unique values of the variable if they number fewer
	than 11. For variables with \eqn{>10} unique values, \code{np}
	equally spaced values in the range are used for plotting if the
	range is not specified.  Variables not specified are set to the default
	adjustment value \code{limits[2]}, i.e. the median for continuous
	variables and a reference category for 	non-continuous ones.
	Later variables define adjustment settings.
	For categorical variables, specify the class labels in quotes when
	specifying variable values.  If the levels of a categorical variable
	are numeric, you may omit the quotes.  For variables not described
	using \code{datadist}, you must specify explicit ranges and
	adjustment settings for predictors  that were in the model.
	If no variables are specified in \dots, predictions will be made by
	separately varying all predictors in the model over their default
	range, holding the other predictors at their adjustment values.
	This has the same effect as specifying \code{name} as a vector
	containing all the predictors.  For \code{rbind}, \dots represents a
	series of results from \code{Predict}.  If you name the results,
	these names will be taken as the values of the new \code{.set.}
	variable added to the concatenated data frames.  See an example below.
  }
  \item{fun}{an optional transformation of the linear predictor}
  \item{type}{
	defaults to providing predictions.  Set to \code{"model.frame"} to
	return a data frame of predictor settings used.  Set to \code{"x"}
	to return the corresponding design matrix constructed from the
	predictor settings.
  }
  \item{np}{
	the number of equally-spaced points computed for continuous
	predictors that vary, i.e., when the specified value is \code{.}
	or \code{NA}
  }
  \item{conf.int}{
	confidence level.  Default is 0.95.  Specify \code{FALSE} to suppress.}
  \item{conf.type}{
	type of confidence interval.  Default is \code{"mean"} which applies
	to all models.  For models containing a residual variance (e.g,
	\code{ols}), you can specify \code{conf.type="individual"} instead,
	to obtain limits on the predicted value for an individual subject.
  }
  \item{adj.zero}{
	Set to \code{TRUE} to adjust all non-plotted variables to 0 (or
	reference cell for categorical variables) and to omit intercept(s)
	from consideration. Default is \code{FALSE}.
  }
  \item{ref.zero}{
	Set to \code{TRUE} to subtract a constant from \eqn{X\beta}{X beta}
	before plotting so that the reference value of the \code{x}-variable
	yields \code{y=0}.  This is done before applying function \code{fun}.
  }
  \item{non.slopes}{
	This is only useful in a multiple intercept model such as the ordinal
	logistic model. There to use to second of three intercepts, for example,
	specify \code{non.slopes=c(0,1,0)}. The default is
	\code{non.slopes=rep(0,k)} if \code{adj.zero=TRUE}, where \code{k}
	is the number of intercepts in the model.  If \code{adj.zero=FALSE},
	the default is \code{(1,0,0,\dots,0)}. 
  }
  \item{time}{
	Specify a single time \code{u} to cause function \code{survest} to
	be invoked to plot the probability of surviving until time \code{u}
	when the fit is from \code{cph} or \code{psm}.
  }
  \item{loglog}{
	Specify \code{loglog=TRUE} to plot \code{log[-log(survival)]}
	instead of survival, when \code{time} is given.
  }
  \item{digits}{
	Controls how ``adjust-to'' values are plotted.  The default is 4
	significant digits.
  }
  \item{name}{
	Instead of specifying the variables to vary in the
	\code{variables} (\dots) list, you can specify one or more variables
	by specifying a vector of character string variable names in the
	\code{name} argument.  Using this mode you cannot specify a list of
	variable values to use; prediction is done as if you had said e.g.
	\code{age} without the equal sign.  Also, interacting factors can
	only be set to their reference values using this notation.
  }
  \item{factors}{
	an alternate way of specifying \dots, mainly for use by
	\code{survplot} or \code{gendata}.  This must be a list with one or
	more values for each variable listed, with \code{NA} values for
	default ranges.
	}
  \item{rename}{
	If you are concatenating predictor sets using \code{rbind} and one
	or more of the variables were renamed for one or more of the sets,
	but these new names represent different versions of the same
	predictors (e.g., using or not using imputation), you can specify a
	named character vector to rename predictors to a central name.  For
	example, specify \code{rename=c(age.imputed='age',
	  corrected.bp='bp')} to rename from old names \code{age.imputed,
	  corrected.bp} to \code{age, bp}.  This happens before
	concatenation of rows.
	}
}
\details{
When there are no intercepts in the fitted model, plot subtracts
adjustment values from each factor while computing variances for
confidence limits. 

Specifying \code{time} will not work for Cox models with time-dependent
covariables.  Use \code{survest} or \code{survfit} for that purpose.
}
\value{
  a data frame containing all model predictors and the computed values
  \code{yhat}, \code{lower}, \code{upper}, the latter two if confidence
  intervals were requested.  The data frame has an additional
  \code{class} \code{"Predict"}.  If \code{name} is specified or no
  predictors are specified in \dots, the resulting data frame has an
  additional variable called \code{.predictor.} specifying which
  predictor is currently being varied.   \code{.predictor.} is handy for
  use as a paneling variable in \code{lattice} or \code{ggplot2} graphics.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{datadist}}, \code{\link{predictrms}},
\code{\link{contrast.rms}}, \code{\link{summary.rms}},  
\code{\link{rms}}, \code{\link{rms.trans}}, \code{\link{survest}},
\code{\link{survplot}}, \code{\link{rmsMisc}},
\code{\link[Hmisc]{transace}}, \code{rbind}
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
label(age)            <- 'Age'      # label is in Hmisc
label(cholesterol)    <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex)            <- 'Sex'
units(cholesterol)    <- 'mg/dl'   # uses units.default in Hmisc
units(blood.pressure) <- 'mmHg'

# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)))
Predict(fit, age, cholesterol, np=4)
Predict(fit, age=seq(20,80,by=10), sex, conf.int=FALSE)
Predict(fit, age=seq(20,80,by=10), sex='male')  # works if datadist not used

ddist$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: ddist$limits["Adjust to","age"] <- 30
fit <- update(fit)   # make new reference value take effect
Predict(fit, age, ref.zero=TRUE, fun=exp)

# Make two curves, and plot the predicted curves as two trellis panels
w <- Predict(fit, age, sex)
require(lattice)
xyplot(yhat ~ age | sex, data=w, type='l')
# To add confidence bands we need to use the Hmisc xYplot function in
# place of xyplot
xYplot(Cbind(yhat,lower,upper) ~ age | sex, data=w, 
       method='filled bands', type='l', col.fill=gray(.95))
# If non-displayed variables were in the model, add a subtitle to show
# their settings using title(sub=paste('Adjusted to',attr(w,'info')$adjust),adj=0)
# Easier: feed w into plot.Predict
\dontrun{
# Predictions form a parametric survival model
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
              rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
t <- -log(runif(n))/h
label(t) <- 'Follow-up Time'
e <- ifelse(t<=cens,1,0)
t <- pmin(t, cens)
units(t) <- "Year"
ddist <- datadist(age, sex)
Srv <- Surv(t,e)

# Fit log-normal survival model and plot median survival time vs. age
f <- psm(Surv(t, e) ~ rcs(age), dist=if(.R.)'lognormal' else 'gaussian')
med <- Quantile(f)       # Creates function to compute quantiles
                         # (median by default)
Predict(f, age, fun=function(x)med(lp=x))
# Note: This works because med() expects the linear predictor (X*beta)
#       as an argument.  Would not work if use 
#       ref.zero=TRUE or adj.zero=TRUE.
# Also, confidence intervals from this method are approximate since
# they don't take into account estimation of scale parameter

# Fit an ols model to log(y) and plot the relationship between x1
# and the predicted mean(y) on the original scale without assuming
# normality of residuals; use the smearing estimator.  Before doing
# that, show confidence intervals for mean and individual log(y)
set.seed(1)
x1 <- runif(300)
x2 <- runif(300)
ddist <- datadist(x1,x2); options(datadist='ddist')
y  <- exp(x1+ x2 - 1 + rnorm(300))
f  <- ols(log(y) ~ pol(x1,2) + x2)
p1 <- Predict(f, x1, conf.type='mean')
p2 <- Predict(f, x1, conf.type='individual')
p  <- rbind(mean=p1, individual=p2)
plot(p, label.curve=FALSE)   # uses superposition
plot(p, ~x1 | .set.)         # 2 panels

r <- resid(f)
smean <- function(yhat)smearingEst(yhat, exp, res, statistic='mean')
formals(smean) <- list(yhat=numeric(0), res=r[!is.na(r)])
#smean$res <- r[!is.na(r)]   # define default res argument to function
Predict(f, x1, fun=smean)
}
options(datadist=NULL)
}
\keyword{models}
