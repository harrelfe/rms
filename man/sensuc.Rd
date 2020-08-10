\name{sensuc}
\alias{sensuc}
\alias{plot.sensuc}
\title{Sensitivity to Unmeasured Covariables}

\description{
Performs an analysis of the sensitivity of a binary treatment (\eqn{X})
effect to an unmeasured binary confounder (\eqn{U}) for a fitted binary
logistic or an unstratified non-time-dependent Cox survival model (the
function works well for the former, not so well for the latter).  This
is done by fitting a sequence of models with separately created \eqn{U}
variables added to the original model.  The sequence of models is formed
by simultaneously varying \eqn{a} and \eqn{b}, where \eqn{a} measures
the association between \eqn{U} and \eqn{X} and \eqn{b} measures the
association between \eqn{U} and \eqn{Y}, where \eqn{Y} is the outcome of
interest.  For Cox models, an approximate solution is used by letting
\eqn{Y} represent some binary classification of the event/censoring time
and the event indicator.  For example, \eqn{Y} could be just be the
event indicator, ignoring time of the event or censoring, or it could be
\eqn{1} if a subject failed before one year and \eqn{0} otherwise.  When
for each combination of \eqn{a} and \eqn{b} the vector of binary values
\eqn{U} is generated, one of two methods is used to constrain the
properties of \eqn{U}.  With either method, the overall prevalence of
\eqn{U} is constrained to be \code{prev.u}.  With the default method
(\code{or.method="x:u y:u"}), \eqn{U} is sampled so that the \eqn{X:U}
odds ratio is \eqn{a} and the \eqn{Y:U} odds ratio is \eqn{b}.  With the
second method, \eqn{U} is sampled according to the model \eqn{logit(U=1
| X, Y) = \alpha + \beta*Y + \gamma*X}, where \eqn{\beta=\log(b)} and
\eqn{\gamma=\log(a)} and \eqn{\alpha} is determined so that the
prevalence of \eqn{U=1} is \code{prev.u}.  This second method results in
the adjusted odds ratio for \eqn{Y:U} given \eqn{X} being \eqn{b}
whereas the default method forces the unconditional (marginal) \eqn{Y:U}
odds ratio to be \eqn{b}.  Rosenbaum uses the default method.

There is a \code{plot} method for plotting objects created by
\code{sensuc}.  Values of \eqn{a} are placed on the x-axis and observed
marginal odds or hazards ratios for \eqn{U} (unadjusted ratios) appear
on the y-axis.  For Cox models, the hazard ratios will not agree exactly
with \eqn{X}:event indicator odds ratios but they sometimes be made
close through judicious choice of the \code{event} function.  The
default plot uses four symbols which differentiate whether for the
\eqn{a,b} combination the effect of \eqn{X} adjusted for \eqn{U} (and
for any other covariables that were in the original model fit) is
positive (usually meaning an effect ratio greater than 1) and
"significant", merely positive, not positive and non significant, or not
positive but significant.  There is also an option to draw the numeric
value of the \eqn{X} effect ratio at the \eqn{a},\eqn{b} combination
along with its \eqn{Z} statistic underneath in smaller letters, and an
option to draw the effect ratio in one of four colors depending on the
significance of the \eqn{Z} statistic.
}
\usage{
# fit <- lrm(formula=y ~ x + other.predictors, x=TRUE, y=TRUE)  #or
# fit <- cph(formula=Surv(event.time,event.indicator) ~ x + other.predictors,
#            x=TRUE, y=TRUE)

sensuc(fit,  
       or.xu=seq(1, 6, by = 0.5), or.u=or.xu, 
       prev.u=0.5, constrain.binary.sample=TRUE, 
       or.method=c("x:u y:u","u|x,y"),
       event=function(y) if(is.matrix(y))y[,ncol(y)] else 1*y)

\method{plot}{sensuc}(x,  ylim=c((1+trunc(min(x$effect.u)-.01))/
                   ifelse(type=='numbers',2,1),
                   1+trunc(max(x$effect.u)-.01)),
     xlab='Odds Ratio for X:U',
     ylab=if(x$type=='lrm')'Odds Ratio for Y:U' else
          'Hazard Ratio for Y:U',
     digits=2, cex.effect=.75, cex.z=.6*cex.effect,
     delta=diff(par('usr')[3:4])/40, 
     type=c('symbols','numbers','colors'),
     pch=c(15,18,5,0), col=c(2,3,1,4), alpha=.05,
     impressive.effect=function(x)x > 1,\dots)
}
\arguments{
\item{fit}{
result of \code{lrm} or \code{cph} with \code{x=TRUE, y=TRUE}.  The
first variable in the right hand side of the model formula must have
been the binary \eqn{X} variable, and it may not interact with other
predictors.
}
\item{x}{
result of \code{sensuc}
}
\item{or.xu}{
vector of possible odds ratios measuring the \eqn{X:U} association.
}
\item{or.u}{
vector of possible odds ratios measuring the \eqn{Y:U} association.
Default is \code{or.xu}.
}
\item{prev.u}{
desired prevalence of \eqn{U=1}.  Default is 0.5, which is usually a
"worst case" for sensitivity analyses.
}
\item{constrain.binary.sample}{
By default, the binary \eqn{U} values are sampled from the appropriate
distributions conditional on \eqn{Y} and \eqn{X} so that the proportions of
\eqn{U=1} in each sample are exactly the desired probabilities, to within
the closeness of \eqn{n\times}probability to an integer.  Specify
\code{constrain.binary.sample=FALSE} to sample from ordinary Bernoulli
distributions, to allow proportions of \eqn{U=1} to reflect sampling fluctuations.
}
\item{or.method}{
see above
}
\item{event}{
a function classifying the response variable into a binary event for the
purposes of constraining the association between \eqn{U} and \eqn{Y}.
For binary logistic models, \code{event} is left at its default value, which
is the identify function, i.e, the original \eqn{Y} values are taken as the
events (no other choice makes any sense here).  For Cox models, the 
default \code{event} function takes the last column of the \code{Surv} object
stored with the fit.  For rare events (high proportion of censored
observations), odds ratios approximate hazard ratios, so the default is OK.  
For other cases, the survival times should be considered (probably in
conjunction with the event indicators), although it may not be possible
to get a high enough hazard ratio between \eqn{U} and \eqn{Y} by sampling \eqn{U} by
temporarily making \eqn{Y} binary.  See the last example which is
for a 2-column \code{Surv} object (first column of
response variable=event time,  second=event indicator).  When
dichotomizing survival time at a given point, it is advantageous to choose
the cutpoint so that not many censored survival times preceed the cutpoint.
Note that in fitting Cox models to examine sensitivity to \eqn{U}, the original
non-dichotomized failure times are used.
}
\item{ylim}{
y-axis limits for \code{plot}
}
\item{xlab}{
x-axis label
}
\item{ylab}{
y-axis label
}
\item{digits}{
number of digits to the right of the decimal point for drawing numbers
on the plot, for
\code{type="numbers"} or \code{type="colors"}.
}
\item{cex.effect}{
character size for drawing effect ratios
}
\item{cex.z}{
character size for drawing \eqn{Z} statistics
}
\item{delta}{
decrement in \eqn{y} value used to draw \eqn{Z} values below effect ratios
}
\item{type}{
specify \code{"symbols"} (the default), \code{"numbers"}, or \code{"colors"} (see above)
}
\item{pch}{
4 plotting characters corresponding to positive and significant
effects for \eqn{X}, positive and non-significant effects, not positive and
not significant, not positive but significant
}
\item{col}{
4 colors as for \code{pch}
}
\item{alpha}{
significance level
}
\item{impressive.effect}{
a function of the odds or hazard ratio for \eqn{X} returning \code{TRUE} for a
positive effect.  By default, a positive effect is taken to mean a
ratio exceeding one.
}
\item{...}{
optional arguments passed to \code{plot}
}}
\value{
\code{sensuc} returns an object of class \code{"sensuc"} with the following elements: \code{OR.xu}
(vector of desired \eqn{X:U} odds ratios or \eqn{a} values), \code{OOR.xu}
(observed marginal \eqn{X:U} odds ratios), \code{OR.u} (desired \eqn{Y:U} odds
ratios or \eqn{b} values), \code{effect.x} (adjusted odds or hazards ratio for
\eqn{X} in a model adjusted for \eqn{U} and all of the other predictors),
\code{effect.u} (unadjusted \eqn{Y:U} odds or hazards ratios), \code{effect.u.adj}
(adjusted \eqn{Y:U} odds or hazards ratios), \eqn{Z} (Z-statistics), \code{prev.u}
(input to \code{sensuc}), \code{cond.prev.u} (matrix with one row per \eqn{a},\eqn{b}
combination, specifying prevalences of \eqn{U} conditional on \eqn{Y} and \eqn{X}
combinations), and \code{type} (\code{"lrm"} or \code{"cph"}).
}
\author{
Frank Harrell\cr
Mark Conaway\cr
Department of Biostatistics\cr
Vanderbilt University School of Medicine\cr
fh@fharrell.com, mconaway@virginia.edu
}
\references{
Rosenbaum, Paul R (1995): Observational Studies.  New York: Springer-Verlag.


Rosenbaum P, Rubin D (1983): Assessing sensitivity to an unobserved binary
covariate in an observational study with binary outcome.  J Roy Statist Soc
B 45:212--218.


Lee WC (2011): Bounding the bias of unmeasured factors with confounding
and effect-modifying potentials.  Stat in Med 30:1007-1017.
}
\seealso{
	\code{\link{lrm}}, \code{\link{cph}}, \code{\link{sample}}
%%	\code{\link[treatSens]{treatSens}}
}
\examples{
set.seed(17)
x <- sample(0:1, 500,TRUE)
y <- sample(0:1, 500,TRUE)
y[1:100] <- x[1:100]  # induce an association between x and y
x2 <- rnorm(500)


f <- lrm(y ~ x + x2, x=TRUE, y=TRUE)


#Note: in absence of U odds ratio for x is exp(2nd coefficient)


g <- sensuc(f, c(1,3))


# Note: If the generated sample of U was typical, the odds ratio for
# x dropped had U been known, where U had an odds ratio
# with x of 3 and an odds ratio with y of 3


plot(g)


# Fit a Cox model and check sensitivity to an unmeasured confounder


# f <- cph(Surv(d.time,death) ~ treatment + pol(age,2)*sex, x=TRUE, y=TRUE)
# sensuc(f, event=function(y) y[,2] & y[,1] < 365.25 )
# Event = failed, with event time before 1 year
# Note: Analysis uses f$y which is a 2-column Surv object
}
\keyword{regression}
\keyword{htest}
\keyword{models}
\keyword{survival}
\concept{model validation}
\concept{sampling}
\concept{logistic regression model}
\concept{sensitivity analysis}

