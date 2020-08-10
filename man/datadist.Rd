\name{datadist}
\alias{datadist}
\alias{print.datadist}
\title{
Distribution Summaries for Predictor Variables
}
\description{
For a given set of variables or a data frame, determines summaries
of variables for effect and plotting ranges, values to adjust to,
and overall ranges
for \code{Predict}, \code{plot.Predict}, \code{ggplot.Predict},
\code{summary.rms}, \code{survplot}, and \code{nomogram.rms}.
If \code{datadist} is called before
a model fit and the resulting object pointed to with \code{options(datadist="name")},
the data characteristics will be stored with the fit by \code{Design()}, so
that later predictions and summaries of the fit will not need to access
the original data used in the fit.  Alternatively, you can specify the
values for each variable in the model when using these 3 functions, or
specify the values of some of them and let the functions look up the
remainder (of say adjustmemt levels) from an object created by \code{datadist}.
The best method is probably to run \code{datadist} once before any models are
fitted, storing the distribution summaries for all potential variables.
Adjustment values are \code{0} for binary variables, the most frequent
category (or optionally the first category level)
for categorical (\code{factor}) variables, the middle level for 
\code{ordered factor} variables, and medians for continuous variables.
See descriptions of \code{q.display} and \code{q.effect} for how display and
effect ranges are chosen for continuous variables.
}
\usage{
datadist(\dots, data, q.display, q.effect=c(0.25, 0.75),
         adjto.cat=c('mode','first'), n.unique=10)

\method{print}{datadist}(x, \dots)
# options(datadist="dd")
# used by summary, plot, survplot, sometimes predict
# For dd substitute the name of the result of datadist
}
\arguments{
\item{...}{
a list of variable names, separated by commas, a single data frame, or
a fit with \code{Design} information.  The first element in this list may
also be an object created by an earlier call to \code{datadist}; then
the later variables are added to this \code{datadist} object.
For a fit object, the variables named
in the fit are retrieved from the active data frame or from the location
pointed to by \code{data=frame number} or \code{data="data frame name"}.
For \code{print}, is ignored.
}
\item{data}{
a data frame or a search position.  If \code{data} is a search position,
it is assumed that a data frame is attached in that position, and all
its variables are used.  If you specify both individual variables in
\code{\dots} and \code{data}, the two sets of variables are combined.  Unless the
first argument is a fit object, \code{data} must be an integer.
}
\item{q.display}{
set of two quantiles for computing the range of continuous variables
to use in displaying regression relationships.  Defaults are
\eqn{q} and \eqn{1-q}, where \eqn{q=10/max(n,200)}, and \eqn{n} is the
number of 
non-missing observations.  Thus for \eqn{n<200}, the .05 and .95 quantiles
are used.  For \eqn{n\geq 200}, the \eqn{10^{th}} smallest and
\eqn{10^{th}} largest values are used.  If you specify \code{q.display},
those quantiles are used whether or not \eqn{n<200}.
}
\item{q.effect}{
set of two quantiles for computing the range of continuous variables
to use in estimating regression effects.  Defaults are c(.25,.75),
which yields inter-quartile-range odds ratios, etc.
}
\item{adjto.cat}{
default is \code{"mode"}, indicating that the modal (most frequent) category
for categorical (factor) variables is the adjust-to setting.
Specify \code{"first"} to use the first level of factor variables as the
adjustment values.  In the case of many levels having the maximum
frequency, the first such level is used for \code{"mode"}.
}
\item{n.unique}{
variables having \code{n.unique} or fewer unique values are considered
to be discrete variables in that their unique values are stored in the
\code{values} list.  This will affect how functions such as
\code{nomogram.Design} determine whether variables are discrete or not.
}
\item{x}{result of \code{datadist}}
}
\value{
a list of class \code{"datadist"} with the following components

\item{limits}{
a \eqn{7 \times k} vector, where \eqn{k} is the number of variables.
The 7 rows correspond to the low value for estimating the effect of
the variable, the value to adjust the variable to when examining
other variables, the high value for effect, low value for displaying
the variable, the high value for displaying it, and the overall lowest
and highest values.
}
\item{values}{
a named list, with one vector of unique values for each numeric
variable having no more than \code{n.unique} unique values
}}
\details{
For categorical variables, the 7 limits are set to character strings
(factors) which correspond to
\code{c(NA,adjto.level,NA,1,k,1,k)}, where \code{k} is the number of levels.
For ordered variables with numeric levels, the limits are set to
\code{c(L,M,H,L,H,L,H)}, where \code{L} is the lowest level, \code{M} is the middle
level, and \code{H} is the highest level.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{rms}}, \code{\link{rms.trans}}, \code{\link[Hmisc]{describe}}, \code{\link{Predict}}, \code{\link{summary.rms}}
}
\examples{
\dontrun{
d <- datadist(data=1)         # use all variables in search pos. 1
d <- datadist(x1, x2, x3)
page(d)                       # if your options(pager) leaves up a pop-up
                              # window, this is a useful guide in analyses
d <- datadist(data=2)         # all variables in search pos. 2
d <- datadist(data=my.data.frame)
d <- datadist(my.data.frame)  # same as previous.  Run for all potential vars.
d <- datadist(x2, x3, data=my.data.frame)   # combine variables
d <- datadist(x2, x3, q.effect=c(.1,.9), q.display=c(0,1))
# uses inter-decile range odds ratios,
# total range of variables for regression function plots
d <- datadist(d, z)           # add a new variable to an existing datadist
options(datadist="d")         #often a good idea, to store info with fit
f <- ols(y ~ x1*x2*x3)


options(datadist=NULL)        #default at start of session
f <- ols(y ~ x1*x2)
d <- datadist(f)              #info not stored in `f'
d$limits["Adjust to","x1"] <- .5   #reset adjustment level to .5
options(datadist="d")


f <- lrm(y ~ x1*x2, data=mydata)
d <- datadist(f, data=mydata)
options(datadist="d")


f <- lrm(y ~ x1*x2)           #datadist not used - specify all values for
summary(f, x1=c(200,500,800), x2=c(1,3,5))         # obtaining predictions
plot(Predict(f, x1=200:800, x2=3))  # or ggplot()


# Change reference value to get a relative odds plot for a logistic model
d$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: d$limits["Adjust to","age"] <- 30
fit <- update(fit)   # make new reference value take effect
plot(Predict(fit, age, ref.zero=TRUE, fun=exp),
     ylab='Age=x:Age=30 Odds Ratio')   # or ggplot()
}
}
\keyword{models}
\keyword{nonparametric}
\keyword{regression}
