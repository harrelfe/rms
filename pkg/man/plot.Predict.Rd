\name{plot.Predict}
\alias{plot.Predict}
\alias{pantext}
\title{Plot Effects of Variables Estimated by a Regression Model Fit}

\description{
  Uses \code{lattice} graphics to plot the effect of one or two predictors
  on the linear predictor or X beta scale, or on some transformation of
  that scale.  The first argument specifies the result of the
  \code{Predict} function.  The predictor is always plotted in its
  original coding.  \code{plot.Predict} uses the
  \code{xYplot} function unless \code{formula} is omitted and the x-axis
  variable is a factor, in which case it reverses the x- and y-axes and
  uses the \code{Dotplot} function.

  If \code{data} is given, a rug plot is drawn showing
  the location/density of data values for the \eqn{x}-axis variable.  If
  there is a \code{groups} (superposition) variable that generated separate
  curves, the data density specific to each class of points is shown.
  This assumes that the second variable was a factor variable.  The rug plots
  are drawn by \code{scat1d}.

  To plot effects instead of estimates (e.g., treatment differences as a
  function of interacting factors) see \code{contrast.rms} and
  \code{summary.rms}.

  \code{pantext} creates a \code{lattice} panel function for including
  text such as that produced by \code{print.anova.rms} inside a panel or
  in a base graphic.
}
\usage{
\method{plot}{Predict}(x, formula, groups=NULL, subset,
     xlim, ylim, xlab, ylab, 
     data=NULL, col.fill=gray(seq(.95, .75, length=5)),
     adj.subtitle, cex.adj, perim=NULL, digits=4, nlevels=3,
     nlines=FALSE, addpanel, scat1d.opts=list(frac=0.025, lwd=0.3), ...)

pantext(object, x, y, cex=.5, adj=0, fontfamily="Courier", lattice=TRUE)

}
\arguments{
  \item{x}{a data frame created by \code{Predict}, or for \code{pantext}
  the x-coordinate for text}
\item{formula}{
  the right hand side of a \code{lattice} formula reference variables in
  data frame \code{x}.  You may not specify \code{formula} if you varied
  multiple predictors separately when calling \code{Predict}.
  Otherwise, when \code{formula} is not given, \code{plot.Predict}
  constructs one from information in \code{x}.
}
\item{groups}{an optional name of one of the variables in \code{x} that
  is to be used as a grouping (superpositioning) variable.  Note that
  \code{groups} does not contain the groups data as is customary in
  \code{lattice}; it is only a single character string specifying the
  name of the grouping variable.}
\item{subset}{a subsetting expression for restricting the rows of
  \code{x} that are used in plotting.  For example, predictions may have
  been requested for males and females but one wants to plot only females.}
\item{xlim}{
This parameter is seldom used, as limits are usually controlled with
\code{Predict}.  One reason to use \code{xlim} is to plot a
\code{factor} variable on the x-axis that was created with the \code{cut2} function
with the \code{levels.mean} option, with \code{val.lev=TRUE} specified to \code{plot.Predict}. 
In this case you may want the axis to
have the range of the original variable values given to \code{cut2} rather
than the range of the means within quantile groups.
}
\item{ylim}{
Range for plotting on response variable axis. Computed by default.
}
\item{xlab}{
Label for \code{x}-axis. Default is one given to \code{asis, rcs}, etc.,
which may have been the \code{"label"} attribute of the variable.
}
\item{ylab}{
Label for \code{y}-axis.  If \code{fun} is not given,
default is \code{"log Odds"} for
\code{lrm}, \code{"log Relative Hazard"} for \code{cph}, name of the response
variable for \code{ols}, \code{TRUE} or \code{log(TRUE)} for \code{psm}, or \code{"X * Beta"} otherwise.
}
\item{data}{a data frame containing the original raw data on which the
  regression model were based, or at least containing the \eqn{x}-axis
  and grouping variable.  If \code{data} is present and contains the
  needed variables, the original data are added to the graph in the form
  of a rug plot using \code{scat1d}.
}
\item{col.fill}{
  a vector of colors used to fill confidence bands for successive
  superposed groups.  Default is inceasingly dark gray scale.
  }
\item{adj.subtitle}{
Set to \code{FALSE} to suppress subtitling the graph with the list of
settings of non-graphed adjustment values.
}
\item{cex.adj}{
\code{cex} parameter for size of adjustment settings in subtitles.  Default is
0.75 times \code{par("cex")}.
}
\item{perim}{
\code{perim} specifies a function having two
arguments.  The first is the vector of values of the first variable that
is about to be plotted on the x-axis.  The second argument is the single
value of the variable representing different curves, for the current
curve being plotted.  The function's returned value must be a logical
vector whose length is the same as that of the first argument, with
values \code{TRUE} if the corresponding point should be plotted for the
current curve, \code{FALSE} otherwise.  See one of the latter examples.
}
\item{digits}{
Controls how numeric variables used for panel labels are formatted. The
default is 4 significant digits.
}
\item{nlevels}{
  when \code{groups} and \code{formula} are not specified, if any panel
  variable has \code{nlevels} or fewer values, that variable is
  converted to a \code{groups} (superpositioning) variable.  Set
  \code{nlevels=0} to prevent this behavior.  For other situations, a
  numeric x-axis variable with \code{nlevels} or fewer unique values
  will cause a dot plot to be drawn instead of an x-y plot.
}
\item{nlines}{If \code{formula} is given, you can set \code{nlines} to
  \code{TRUE} to convert the x-axis variable to a factor and then to an
  integer.  Points are plotted at integer values on the x-axis but
  labeled with category levels.  Points are connected by lines.}
\item{addpanel}{an additional panel function to call along with panel
  functions used for \code{xYplot} and \code{Dotplot} displays}
\item{scat1d.opts}{a list containing named elements that specifies
  parameters to \code{\link{scat1d}} when \code{data} is given.  The
  \code{col} parameter is usually derived from other plotting
  information and not specified by the user.} 
\item{\dots}{
  extra arguments to pass to \code{xYplot} or \code{Dotplot}.  Some
  useful ones are \code{label.curves} and \code{abline}.
  Set \code{label.curves} to \code{FALSE} to suppress labeling of
  separate curves. Default is \code{TRUE}, which
  causes \code{labcurve} to be invoked to place labels at positions where the
  curves are most separated, labeling each curve with the full curve label.
  Set \code{label.curves} to a \code{list} to specify options to
  \code{labcurve}, e.g., \code{label.curves=} \code{list(method="arrow",
	cex=.8)}. 
  These option names may be abbreviated in the usual way arguments
  are abbreviated.  Use for example \code{label.curves=list(keys=letters[1:5])}
  to draw single lower case letters on 5 curves where they are most
  separated, and automatically position a legend
  in the most empty part of the plot.  The \code{col}, \code{lty}, and
  \code{lwd} parameters are passed automatically to \code{labcurve}
  although they may be overridden here.
}
\item{object}{an object having a \code{print} method}
\item{y}{y-coordinate for placing text in a \code{lattice} panel
or on a base graphics plot}
\item{cex}{character expansion size for \code{pantext}}
\item{adj}{text justification.  Default is left justified.}
\item{fontfamily}{
font family for \code{pantext}.  Default is \code{"Courier"} which
will line up columns of a table.
}
\item{lattice}{set to \code{FALSE} to use \code{text} instead of
  \code{ltext} in the function generated by \code{pantext}, to use base
  graphics}
}
\value{
  a \code{lattice} object ready to \code{print} for rendering.
}
\details{
When a \code{groups} (superpositioning) variable was used, you can issue
the command \code{Key(\dots)} after printing the result of
\code{plot.Predict}, to draw a key for the groups.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
  Fox J, Hong J (2009): Effect displays in R for multinomial and
  proportional-odds logit models: Extensions to the effects package.  J
  Stat Software 32 No. 1.
  }
\seealso{
  \code{\link{Predict}}, \code{\link{rbind.Predict}},
  \code{\link{datadist}}, \code{\link{predictrms}}, \code{\link{anova.rms}},
  \code{\link{contrast.rms}}, \code{\link{summary.rms}},
  \code{\link{rms}}, \code{\link{rmsMisc}}, 
  \code{\link[Hmisc]{labcurve}}, \code{\link[Hmisc]{scat1d}},
  \code{\link[Hmisc]{xYplot}}, \code{\link[Hmisc]{Overview}}
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

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)
plot(Predict(fit))       # Plot effects of all 4 predictors
plot(Predict(fit), data=llist(blood.pressure,age))
                         # rug plot for two of the predictors

p <- Predict(fit, name=c('age','cholesterol'))   # Make 2 plots
plot(p)

p <- Predict(fit, age=seq(20,80,length=100), sex, conf.int=FALSE)
                         # Plot relationship between age and log
                         # odds, separate curve for each sex,
plot(p)                  # no C.I.

p <- Predict(fit, age, sex)
plot(p, label.curves=FALSE, data=llist(age,sex))
                         # use label.curves=list(keys=c('a','b'))'
                         # to use 1-letter abbreviations
                         # data= allows rug plots (1-dimensional scatterplots)
                         # on each sex's curve, with sex-
                         # specific density of age
                         # If data were in data frame could have used that
p <- Predict(fit, age=seq(20,80,length=100), sex='male', fun=plogis)
                         # works if datadist not used
plot(p, ylab=expression(hat(P)))
                         # plot predicted probability in place of log odds

per <- function(x, y) x >= 30
plot(p, perim=per)       # suppress output for age < 30 but leave scale alone

# Take charge of the plot setup by specifying a lattice formula
p <- Predict(fit, age, blood.pressure=c(120,140,160),
             cholesterol=c(180,200,215), sex)
plot(p, ~ age | blood.pressure*cholesterol, subset=sex=='male')
plot(p, ~ age | cholesterol*blood.pressure, subset=sex=='female')
plot(p, ~ blood.pressure|cholesterol*round(age,-1), subset=sex=='male')
plot(p)

# Plot the age effect as an odds ratio
# comparing the age shown on the x-axis to age=30 years

ddist$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: ddist$limits["Adjust to","age"] <- 30
fit <- update(fit)   # make new reference value take effect
p <- Predict(fit, age, ref.zero=TRUE, fun=exp)
plot(p, ylab='Age=x:Age=30 Odds Ratio',
     abline=list(list(h=1, lty=2, col=2), list(v=30, lty=2, col=2)))


# Plots for a parametric survival model
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
f <- psm(Surv(t, e) ~ rcs(age), dist='lognormal')
med <- Quantile(f)       # Creates function to compute quantiles
                         # (median by default)
p <- Predict(f, age, fun=function(x) med(lp=x))
plot(p, ylab="Median Survival Time")
# Note: confidence intervals from this method are approximate since
# they don't take into account estimation of scale parameter


# Fit an ols model to log(y) and plot the relationship between x1
# and the predicted mean(y) on the original scale without assuming
# normality of residuals; use the smearing estimator
# See help file for rbind.Predict for a method of showing two
# types of confidence intervals simultaneously.
set.seed(1)
x1 <- runif(300)
x2 <- runif(300)
ddist <- datadist(x1,x2)
y  <- exp(x1+x2-1+rnorm(300))
f <- ols(log(y) ~ pol(x1,2)+x2)
r <- resid(f)
smean <- function(yhat)smearingEst(yhat, exp, res, statistic='mean')
formals(smean) <- list(yhat=numeric(0), res=r[!is.na(r)])
#smean$res <- r[!is.na(r)]   # define default res argument to function
plot(Predict(f, x1, fun=smean), ylab='Predicted Mean on y-scale')

# Make an 'interaction plot', forcing the x-axis variable to be
# plotted at integer values but labeled with category levels
n <- 100
set.seed(1)
gender <- c(rep('male', n), rep('female',n))
m <- sample(c('a','b'), 2*n, TRUE)
d <-  datadist(gender, m); options(datadist='d')
anxiety <- runif(2*n) + .2*(gender=='female') + .4*(gender=='female' & m=='b')
tapply(anxiety, llist(gender,m), mean)
f <- ols(anxiety ~ gender*m)
p <- Predict(f, gender, m)
plot(p)     # horizontal dot chart; usually preferred for categorical predictors
Key(.5, .5)
plot(p, ~gender, groups='m', nlines=TRUE)
plot(p, ~m, groups='gender', nlines=TRUE)
plot(p, ~gender|m, nlines=TRUE)

options(datadist=NULL)

\dontrun{
# Example in which separate curves are shown for 4 income values
# For each curve the estimated percentage of voters voting for
# the democratic party is plotted against the percent of voters
# who graduated from college.  Data are county-level percents.

incomes <- seq(22900, 32800, length=4)  
# equally spaced to outer quintiles
p <- Predict(f, college, income=incomes, conf.int=FALSE)
plot(p, xlim=c(0,35), ylim=c(30,55))

# Erase end portions of each curve where there are fewer than 10 counties having
# percent of college graduates to the left of the x-coordinate being plotted,
# for the subset of counties having median family income with 1650
# of the target income for the curve

show.pts <- function(college.pts, income.pt) {
  s <- abs(income - income.pt) < 1650  #assumes income known to top frame
  x <- college[s]
  x <- sort(x[!is.na(x)])
  n <- length(x)
  low <- x[10]; high <- x[n-9]
  college.pts >= low & college.pts <= high
}

plot(p, xlim=c(0,35), ylim=c(30,55), perim=show.pts)
}
}
\keyword{models}
\keyword{hplot}
\keyword{htest}
