\name{ggplot.Predict}
\alias{ggplot.Predict}
\title{Plot Effects of Variables Estimated by a Regression Model Fit
	Using ggplot2}

\description{
  Uses \code{ggplot2} graphics to plot the effect of one or two predictors
  on the linear predictor or X beta scale, or on some transformation of
  that scale.  The first argument specifies the result of the
  \code{Predict} function.  The predictor is always plotted in its
  original coding.

  If \code{rdata} is given, a spike histogram is drawn showing
  the location/density of data values for the \eqn{x}-axis variable.  If
  there is a \code{groups} (superposition) variable that generated separate
  curves, the data density specific to each class of points is shown.
  This assumes that the second variable was a factor variable.  The histograms
  are drawn by \code{histSpikeg}.

  To plot effects instead of estimates (e.g., treatment differences as a
  function of interacting factors) see \code{contrast.rms} and
  \code{summary.rms}.
}
\usage{
\method{ggplot}{Predict}(data, mapping, formula=NULL, groups=NULL,
     aestype=c('color', 'linetype'),
     conf=c('fill', 'lines'),
     varypred=FALSE, sepdiscrete=c('no', 'list', 'vertical', 'horizontal'),
     subset, xlim., ylim., xlab, ylab, 
     colorscale=function(...) scale_color_manual(...,
       values=c("#000000", "#E69F00", "#56B4E9",
                "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")),
     colfill='black',
     rdata=NULL, anova=NULL, pval=FALSE, size.anova=4,
     adj.subtitle, size.adj=2.5, perim=NULL, nlevels=3,
     flipxdiscrete=TRUE,
     legend.position='right', legend.label=NULL,
     vnames=c('labels','names'), abbrev=FALSE, minlength=6,
     layout=NULL, addlayer,
     histSpike.opts=list(frac=function(f) 0.01 + 
         0.02 * sqrt(f - 1)/sqrt(max(f, 2) - 1), side=1, nint=100),
     type=NULL, ggexpr=FALSE, ..., environment)
}
\arguments{
  \item{data}{a data frame created by \code{Predict}}
	\item{mapping}{kept because of \code{ggplot} generic setup.  If
  specified it will be assumed to be \code{formula}.}
	\item{formula}{
		a \code{ggplot} faceting formula of the form \code{vertical variables
		~ horizontal variables}, with variables separated by \code{*} if
	there is more than one variable on a side.  If omitted, the formula
	will be built using assumptions on the list of variables that varied
	in the \code{Predict} call.  When plotting multiple panels (for
  separate predictors), \code{formula} may be specified but by default
  no formula is constructed.
  }
\item{groups}{an optional character string containing the
	name of one of the variables in \code{data} that
  is to be used as a grouping (superpositioning) variable.
  Set \code{groups=FALSE} to suppress superpositioning.  By default, the
  second varying variable is used for superpositioning groups.  You can
  also specify a length 2 string vector of variable names specifying two
  dimensions of superpositioning, identified by different aesthetics
  corresponding to the \code{aestype} argument. When plotting effects
  of more than one predictor, \code{groups}
  is a character string that specifies a single variable name in
  \code{data} that can be used to form panels.  Only applies if using
  \code{rbind} to combine several \code{Predict} results. }
\item{aestype}{a string vector of aesthetic names corresponding to
  variables in the \code{groups} vector.  Default is to use, in order,
  \code{color}, and \code{linetype}.  Other permissible values are
	\code{size}, \code{shape}.}
\item{conf}{specify \code{conf="line"} to show confidence bands with
  lines instead of filled ribbons, the default}
\item{varypred}{set to \code{TRUE} if \code{data} is the result of
  passing multiple \code{Predict} results, that represent different
  predictors, to \code{rbind.Predict}.  This will cause the \code{.set.}
  variable created by \code{rbind} to be copied to the
  \code{.predictor.} variable.}
\item{sepdiscrete}{set to something other than \code{"no"} to create
  separate graphics for continuous and discrete predictors.  For
  discrete predictors, horizontal dot charts are produced.  This allows
  use of the \code{ggplot2} \code{facet_wrap} function to make better
  use of space.  If \code{sepdiscrete="list"}, a list of two \code{grid}
  graphics objects is returned if both types of predictors are present
  (otherwise one object for the type that existed in the model).  Set
  \code{sepdiscrete="vertical"} to put the two types of plots into one
  graphical object with continuous predictors on top and given a
  fraction of space relative to the number of continuous vs. number of
  discrete variables.  Set \code{sepdiscrete="horizontal"} to get a
  horizontal arrangements with continuous variables on the left.}
\item{subset}{a subsetting expression for restricting the rows of
  \code{data} that are used in plotting.  For example, predictions may have
  been requested for males and females but one wants to plot only females.}
\item{xlim.}{
This parameter is seldom used, as limits are usually controlled with
\code{Predict}.  Usually given as its legal abbreviation \code{xlim}.
  One reason to use \code{xlim} is to plot a \code{factor} variable on
  the x-axis that was created with the \code{cut2} function with the
  \code{levels.mean} option, with \code{val.lev=TRUE} specified to
  \code{plot.Predict}.  In this case you may want the axis to have the
  range of the original variable values given to \code{cut2} rather than
  the range of the means within quantile groups. 
}
\item{ylim.}{
Range for plotting on response variable axis. Computed by default.
  Usually specified using its legal definition \code{ylim}.
}
\item{xlab}{
Label for \code{x}-axis. Default is one given to \code{asis, rcs}, etc.,
which may have been the \code{"label"} attribute of the variable.
}
\item{ylab}{
Label for \code{y}-axis.  If \code{fun} is not given,
default is \code{"log Odds"} for
\code{lrm}, \code{"log Relative Hazard"} for \code{cph}, name of the response
variable for \code{ols}, \code{TRUE} or \code{log(TRUE)} for \code{psm},
  or \code{"X * Beta"} otherwise.  Specify \code{ylab=NULL} to omit
  \code{y}-axis labels.
}
\item{colorscale}{a \code{ggplot2} discrete scale function,
  e.g. \code{function(...) scale_color_brewer(..., palette='Set1',
  type='qual')}.  The default is the colorblind-friendly palette
  including black in \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)}.
}
\item{colfill}{a single character string or number specifying the fill color
  to use for \code{geom_ribbon} for shaded confidence bands.  Alpha
  transparency of 0.2 is applied to any color specified.}
\item{rdata}{a data frame containing the original raw data on which the
  regression model were based, or at least containing the \eqn{x}-axis
  and grouping variable.  If \code{rdata} is present and contains the
  needed variables, the original data are added to the graph in the form
  of a spike histogram using \code{histSpikeg} in the Hmisc package.
}
\item{anova}{an object returned by \code{\link{anova.rms}}.  If
	\code{anova} is specified, the overall test of association for
	predictor plotted is added as text to each panel, located at the spot
	at which the panel is most empty unless there is significant empty
	space at the top or bottom of the panel; these areas are given preference.}
\item{pval}{specify \code{pval=TRUE} for \code{anova} to include not
	only the test statistic but also the P-value}
\item{size.anova}{character size for the test statistic printed on the
	panel, mm}
\item{adj.subtitle}{
Set to \code{FALSE} to suppress subtitling the graph with the list of
settings of non-graphed adjustment values.
}
\item{size.adj}{Size of adjustment settings in subtitles in mm.  Default is 2.5.}
\item{perim}{
\code{perim} specifies a function having two
arguments.  The first is the vector of values of the first variable that
is about to be plotted on the x-axis.  The second argument is the single
value of the variable representing different curves, for the current
curve being plotted.  The function's returned value must be a logical
vector whose length is the same as that of the first argument, with
values \code{TRUE} if the corresponding point should be plotted for the
current curve, \code{FALSE} otherwise.  See one of the latter examples.
\code{perim} only applies if predictors were specified to \code{Predict}.
}
\item{nlevels}{
  when \code{groups} and \code{formula} are not specified, if any panel
  variable has \code{nlevels} or fewer values, that variable is
  converted to a \code{groups} (superpositioning) variable.  Set
  \code{nlevels=0} to prevent this behavior.  For other situations, a
  non-numeric x-axis variable with \code{nlevels} or fewer unique values
  will cause a horizontal dot plot to be drawn instead of an x-y plot
  unless \code{flipxdiscrete=FALSE}.
}
\item{flipxdiscrete}{see \code{nlevels}}
\item{legend.position}{\code{"right"} (the default for single-panel
  plots), \code{"left"}, \code{"bottom"}, \code{"top"}, a two-element
  numeric vector, or \code{"none"} to suppress.  For multi-panel plots
  the default is \code{"top"}, and a legend only appears for the first
  (top left) panel.}
\item{legend.label}{if omitted, group variable labels will be used for
  label the legend.  Specify \code{legend.label=FALSE} to suppress using
  a legend name, or a character string or expression to specify the
  label.  Can be a vector is there is more than one grouping variable.}
\item{vnames}{applies to the case where multiple plots are produced
  separately by predictor.  Set to \code{'names'} to use variable names
  instead of labels for these small plots.}
\item{abbrev}{set to true to abbreviate levels of predictors that are
  categorical to a minimum length of \code{minlength}}
\item{minlength}{see \code{abbrev}}
\item{layout}{for multi-panel plots a 2-vector specifying the number of
  rows and number of columns.  If omitted will be computed from the
  number of panels to make as square as possible.}
\item{addlayer}{a \code{ggplot2} expression consisting of one or more
  layers to add to the current plot}
\item{histSpike.opts}{a list containing named elements that specifies
  parameters to \code{\link{histSpikeg}} when \code{rdata} is given.  The
  \code{col} parameter is usually derived from other plotting
  information and not specified by the user.}
\item{type}{a value (\code{"l","p","b"}) to override default choices
  related to showing or connecting points.  Especially  useful for
  discrete x coordinate variables.}
\item{ggexpr}{set to \code{TRUE} to have the function return the
  character string(s) constructed to invoke \code{ggplot} without
  executing the commands}
\item{\dots}{ignored}
\item{environment}{ignored; used to satisfy rules because of the generic ggplot}
}
\value{an object of class \code{"ggplot2"} ready for printing.  For the
  case where predictors were not specified to \code{Predict}, 
  \code{sepdiscrete=TRUE}, and there were both continuous and discrete
  predictors in the model, a list of two graphics objects is returned.}
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
\note{If plotting the effects of all predictors you can reorder the
  panels using for example \code{p <- Predict(fit); p$.predictor. <-
	factor(p$.predictor., v)} where \code{v} is a vector of predictor
  names specified in the desired order.
  }
\seealso{
  \code{\link{Predict}}, \code{\link{rbind.Predict}},
  \code{\link{datadist}}, \code{\link{predictrms}}, \code{\link{anova.rms}},
  \code{\link{contrast.rms}}, \code{\link{summary.rms}},
  \code{\link{rms}}, \code{\link{rmsMisc}}, \code{\link{plot.Predict}},
  \code{\link[Hmisc]{labcurve}}, \code{\link[Hmisc]{histSpikeg}},
  \code{\link[ggplot2]{ggplot}}, \code{\link[Hmisc]{Overview}}
}
\examples{
n <- 350     # define sample size
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
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male')) +
  .01 * (blood.pressure - 120)
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)
an <- anova(fit)
# Plot effects in two vertical sub-panels with continuous predictors on top
# ggplot(Predict(fit), sepdiscrete='vertical')
# Plot effects of all 4 predictors with test statistics from anova, and P
ggplot(Predict(fit), anova=an, pval=TRUE)
# ggplot(Predict(fit), rdata=llist(blood.pressure, age))
# spike histogram plot for two of the predictors

# p <- Predict(fit, name=c('age','cholesterol'))   # Make 2 plots
# ggplot(p)

# p <- Predict(fit, age=seq(20,80,length=100), sex, conf.int=FALSE)
#                        # Plot relationship between age and log
                         # odds, separate curve for each sex,
# ggplot(p, subset=sex=='female' | age > 30)
# No confidence interval, suppress estimates for males <= 30

# p <- Predict(fit, age, sex)
# ggplot(p, rdata=llist(age,sex))
                         # rdata= allows rug plots (1-dimensional scatterplots)
                         # on each sex's curve, with sex-
                         # specific density of age
                         # If data were in data frame could have used that
# p <- Predict(fit, age=seq(20,80,length=100), sex='male', fun=plogis)
                         # works if datadist not used
# ggplot(p, ylab=expression(hat(P)))
                         # plot predicted probability in place of log odds
# per <- function(x, y) x >= 30
# ggplot(p, perim=per)       # suppress output for age < 30 but leave scale alone

# Do ggplot2 faceting a few different ways
p <- Predict(fit, age, sex, blood.pressure=c(120,140,160),
             cholesterol=c(180,200,215))
# ggplot(p)
ggplot(p, cholesterol ~ blood.pressure)
# ggplot(p, ~ cholesterol + blood.pressure)
# color for sex, line type for blood.pressure:
ggplot(p, groups=c('sex', 'blood.pressure'))
# Add legend.position='top' to allow wider plot
# Map blood.pressure to line thickness instead of line type:
# ggplot(p, groups=c('sex', 'blood.pressure'), aestype=c('color', 'size'))

# Plot the age effect as an odds ratio
# comparing the age shown on the x-axis to age=30 years

# ddist$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: ddist$limits["Adjust to","age"] <- 30
# fit <- update(fit)   # make new reference value take effect
# p <- Predict(fit, age, ref.zero=TRUE, fun=exp)
# ggplot(p, ylab='Age=x:Age=30 Odds Ratio',
#        addlayer=geom_hline(yintercept=1, col=gray(.8)) +
#                 geom_vline(xintercept=30, col=gray(.8)) +
#                 scale_y_continuous(trans='log',
#                       breaks=c(.5, 1, 2, 4, 8))))

# Compute predictions for three predictors, with superpositioning or
# conditioning on sex, combined into one graph

p1 <- Predict(fit, age, sex)
p2 <- Predict(fit, cholesterol, sex)
p3 <- Predict(fit, blood.pressure, sex)
p <- rbind(age=p1, cholesterol=p2, blood.pressure=p3)
ggplot(p, groups='sex', varypred=TRUE, adj.subtitle=FALSE)
# ggplot(p, groups='sex', varypred=TRUE, adj.subtitle=FALSE, sepdiscrete='vert')

\dontrun{
# For males at the median blood pressure and cholesterol, plot 3 types
# of confidence intervals for the probability on one plot, for varying age
ages <- seq(20, 80, length=100)
p1 <- Predict(fit, age=ages, sex='male', fun=plogis)  # standard pointwise
p2 <- Predict(fit, age=ages, sex='male', fun=plogis,
              conf.type='simultaneous')               # simultaneous
p3 <- Predict(fit, age=c(60,65,70), sex='male', fun=plogis,
              conf.type='simultaneous')               # simultaneous 3 pts
# The previous only adjusts for a multiplicity of 3 points instead of 100
f <- update(fit, x=TRUE, y=TRUE)
g <- bootcov(f, B=500, coef.reps=TRUE)
p4 <- Predict(g, age=ages, sex='male', fun=plogis)    # bootstrap percentile
p <- rbind(Pointwise=p1, 'Simultaneous 100 ages'=p2,
           'Simultaneous     3 ages'=p3, 'Bootstrap nonparametric'=p4)
# as.data.frame so will call built-in ggplot
ggplot(as.data.frame(p), aes(x=age, y=yhat)) + geom_line() +
 geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2, linetype=0)+
 facet_wrap(~ .set., ncol=2)

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
f <- psm(Srv ~ rcs(age), dist='lognormal')
med <- Quantile(f)       # Creates function to compute quantiles
                         # (median by default)
p <- Predict(f, age, fun=function(x) med(lp=x))
ggplot(p, ylab="Median Survival Time")
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
ggplot(Predict(f, x1, fun=smean), ylab='Predicted Mean on y-scale')
}

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
# ggplot(p)     # horizontal dot chart; usually preferred for categorical predictors
# ggplot(p, flipxdiscrete=FALSE)  # back to vertical
ggplot(p, groups='gender')
ggplot(p, ~ m, groups=FALSE, flipxdiscrete=FALSE)

options(datadist=NULL)

\dontrun{
# Example in which separate curves are shown for 4 income values
# For each curve the estimated percentage of voters voting for
# the democratic party is plotted against the percent of voters
# who graduated from college.  Data are county-level percents.

incomes <- seq(22900, 32800, length=4)  
# equally spaced to outer quintiles
p <- Predict(f, college, income=incomes, conf.int=FALSE)
ggplot(p, xlim=c(0,35), ylim=c(30,55))

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

ggplot(p, xlim=c(0,35), ylim=c(30,55), perim=show.pts)

# Rename variables for better plotting of a long list of predictors
f <- ...
p <- Predict(f)
re <- c(trt='treatment', diabet='diabetes', sbp='systolic blood pressure')

for(n in names(re)) {
  names(p)[names(p)==n] <- re[n]
  p$.predictor.[p$.predictor.==n] <- re[n]
  }
ggplot(p)
}
}
\keyword{models}
\keyword{hplot}
\keyword{htest}
