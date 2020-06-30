\name{nomogram}
\alias{nomogram}
\alias{print.nomogram}
\alias{plot.nomogram}
\alias{legend.nomabbrev}
\title{Draw a Nomogram Representing a Regression Fit}
\description{
  Draws a partial nomogram that can be used to manually obtain predicted
  values from a regression model that was fitted with \code{rms}.
  The nomogram does not have lines representing sums, but it has a reference
  line for reading scoring points (default range 0--100).  Once the reader
  manually totals the points, the predicted values can be read at the bottom.
  Non-monotonic transformations of continuous variables are handled (scales
  wrap around), as
  are transformations which have flat sections (tick marks are labeled
  with ranges).  If interactions are in the model, one variable
  is picked as the \dQuote{axis variable}, and separate axes are constructed for
  each level of the interacting factors (preference is given automatically
  to using any discrete factors to construct separate axes) and
  levels of factors which are indirectly related to interacting
  factors (see DETAILS).  Thus the nomogram is designed so that only
  one axis is actually read for each variable, since the variable
  combinations are disjoint.  For
  categorical interacting factors, the default is to construct axes for
  all levels.
  The user may specify
  coordinates of each predictor to label on its axis, or use default values.
  If a factor interacts with other factors, settings for one or more of
  the interacting factors may be specified separately (this is mandatory
  for continuous variables).  Optional confidence intervals will be
  drawn for individual scores as well as for the linear predictor.
  If more than one confidence level is chosen, multiple levels may be
  displayed using different colors or gray scales.  Functions of the
  linear predictors may be added to the nomogram.

	The \code{\link{datadist}} object that was in effect when the model
	was fit is used to specify the limits of the axis for continuous
	predictors when the user does not specify tick mark locations in the
	\code{nomogram} call.

  \code{print.nomogram} prints axis information stored in an object returned
  by \code{nomogram}.  This is useful in producing tables of point assignments
  by levels of predictors.  It also prints how many linear predictor
  units there are per point and the number of points per unit change in
  the linear predictor.

  \code{legend.nomabbrev} draws legends describing abbreviations used for
  labeling tick marks for levels of categorical predictors.
}
\usage{
nomogram(fit, ..., adj.to, lp=TRUE, lp.at=NULL,
         fun=NULL, fun.at=NULL, fun.lp.at=NULL, funlabel="Predicted Value",
         interact=NULL, kint=NULL,  conf.int=FALSE, 
         conf.lp=c("representative", "all", "none"),
         est.all=TRUE, posterior.summary=c('mean', 'median', 'mode'),
         abbrev=FALSE, minlength=4, maxscale=100, nint=10, 
         vnames=c("labels","names"),
         varname.label=TRUE, varname.label.sep="=",
         omit=NULL, verbose=FALSE)

\method{print}{nomogram}(x, dec=0, \dots)

\method{plot}{nomogram}(x, lplabel="Linear Predictor", fun.side,
 col.conf=c(1, 0.3),
 conf.space=c(.08,.2), label.every=1, force.label=FALSE, 
 xfrac=.35, cex.axis=.85, cex.var=1, col.grid=NULL,
 varname.label=TRUE, varname.label.sep="=", ia.space=.7, 
 tck=NA, tcl=-0.25, lmgp=.4, naxes,
 points.label='Points', total.points.label='Total Points',
 total.sep.page=FALSE, total.fun, cap.labels=FALSE, ...)

legend.nomabbrev(object, which, x, y, ncol=3, \dots)
}
\arguments{
  \item{fit}{
    a regression model fit that was created with \code{rms}, and
	(usually) with \code{options(datadist = "object.name")} in effect. 
  }
  \item{\dots}{
    settings of variables to use in constructing axes.  If \code{datadist}
    was in effect, the default is to use \code{pretty(total range, nint)}
    for continuous variables, and the class levels for discrete ones. 
    For \code{legend.nomabbrev}, \code{\dots} specifies optional
	parameters to pass 
    to \code{legend}.  Common ones are \code{bty = "n"} to suppress drawing the
    box.  You may want to specify a non-proportionally spaced font
    (e.g., courier) number if abbreviations are more than one letter long.
    This will make the abbreviation definitions line up (e.g., specify
    \code{font = 2}, the default for courier).  Ignored for \code{print}
	and \code{plot}.
  }
  \item{adj.to}{
    If you didn't define \code{datadist} for all predictors, you will have to
    define adjustment settings for the undefined ones, e.g.
    \code{adj.to= list(age = 50, sex = "female")}.
  }
  \item{lp}{
    Set to \code{FALSE} to suppress creation of an axis for scoring
    \eqn{X\beta}{X beta}.
  }
  \item{lp.at}{
    If \code{lp=TRUE}, \code{lp.at} may specify a vector of settings of
    \eqn{X\beta}{X beta}.
    Default is to use \code{pretty(range of linear predictors, nint)}.
  }
  \item{fun}{
    an optional function to transform the linear predictors, and to plot
    on another axis.  If more than one transformation is plotted, put
    them in a list, e.g. \code{list(function(x) x/2, function(x) 2*x)}.
    Any function values equal to \code{NA} will be ignored.  
  }
  \item{fun.at}{
    function values to label on axis.  Default \code{fun} evaluated
    at \code{lp.at}.   If more than one \code{fun} was specified, using a vector
    for \code{fun.at} will cause all functions to be evaluated at the same
    argument values.  To use different values, specify a list of vectors for
    \code{fun.at}, with elements corresponding to the different functions
    (lists of vectors also applies to \code{fun.lp.at} and \code{fun.side}).
  }
  \item{fun.lp.at}{
    If you want to
    evaluate one of the functions at a different set of linear predictor
    values than may have been used in constructing the linear predictor axis,
    specify a vector or list of vectors 
    of linear predictor values at which to evaluate the function.  This is
    especially useful for discrete functions.  The presence of this attribute
    also does away with the need for \code{nomogram} to compute numerical approximations of 
    the inverse of the function.  It also allows the user-supplied function
    to return \code{factor} objects, which is useful when e.g. a single tick
    mark position actually represents a range.
    If the \code{fun.lp.at} parameter is present, the \code{fun.at}
    vector for that function is ignored.
  }
  \item{funlabel}{
    label for \code{fun} axis.  If more than one function was given but
    funlabel is of length one, it will be duplicated as needed.  If \code{fun} is
    a list of functions for which you specified names (see the final example
    below), these names will be used as labels.
  }
  \item{interact}{
    When a continuous variable interacts with a discrete one, axes are
    constructed so that the continuous variable moves within the axis, and
    separate axes represent levels of interacting factors.  For interactions
    between two continuous variables, all but the axis variable must have
    discrete levels defined in \code{interact}.  
    For discrete interacting factors, you may specify levels to use in
    constructing the multiple axes.  For continuous interacting factors,
    you must do this.  Examples: \code{interact = list(age = seq(10,70,by=10),
      treat = c("A","B","D"))}.
  }
  \item{kint}{
    for models such as the ordinal models with multiple intercepts,
    specifies which one to use in evaluating the linear predictor.
		Default is to use \code{fit$interceptRef} if it exists, or 1.
  }
  \item{conf.int}{
    confidence levels to display for each scoring.  Default is \code{FALSE} to display
    no confidence limits.  Setting \code{conf.int} to \code{TRUE} is the same as
    setting it to \code{c(0.7, 0.9)},
    with the line segment between the 0.7 and 0.9 levels shaded using
    gray scale.  
  }
  \item{conf.lp}{
    default is \code{"representative"} to group all linear predictors evaluated
    into deciles, and to show, for the linear predictor confidence intervals,
    only the mean linear predictor within the deciles along with the median
    standard error within the deciles.  Set \code{conf.lp = "none"} to suppress
    confidence limits for the linear predictors, and to \code{"all"} to show
    all confidence limits.
  }
  \item{est.all}{
    To plot axes for only the subset of variables named in \code{\dots}, set
    \code{est.all = FALSE}.  Note: This option only works when zero has a special
    meaning for the variables that are omitted from the graph.
  }
	\item{posterior.summary}{when operating on a Bayesian model such as a
		result of \code{blrm} specifies whether to use posterior mean
		(default) vs. posterior mode/median of parameter values in constructing
		the nomogram}
  \item{abbrev}{
    Set to \code{TRUE} to use the \code{\link{abbreviate}} function to abbreviate levels of
    categorical factors, both for labeling tick marks and for axis titles.
    If you only want to abbreviate certain predictor variables, set \code{abbrev}
    to a vector of character strings containing their names.
  }
  \item{minlength}{
    applies if \code{abbrev = TRUE}.  Is the minimum abbreviation length passed to the
    \code{\link{abbreviate}} function.  If you set \code{minlength = 1}, the letters of the
    alphabet are used to label tick marks for categorical predictors, and
    all letters are drawn no matter how close together they are.  For
    labeling axes (interaction settings), \code{minlength = 1} causes
    \code{minlength = 4} to be used.
  }
  \item{maxscale}{
    default maximum point score is 100
  }
  \item{nint}{
    number of intervals to label for axes representing continuous variables.
    See \code{\link{pretty}}.
  }
  \item{vnames}{
    By default, variable labels are used to label axes.  Set
    \code{vnames = "names"}
    to instead use variable names.
  }
  \item{omit}{
    vector of character strings containing names of variables for which to
    suppress drawing axes.  Default is to show all variables.
  }
  \item{verbose}{
    set to \code{TRUE} to get printed output detailing how tick marks are chosen
    and labeled for function axes.  This is useful in seeing how certain
    linear predictor values cannot be solved for using inverse linear
    interpolation on the (requested linear predictor values, function values at 
    these lp values).  When this happens you will see \code{NA}s in the verbose
    output, and the corresponding tick marks will not appear in the nomogram.
  }
  \item{x}{an object created by \code{nomogram}, or the x coordinate for
	a legend}
  \item{dec}{
    number of digits to the right of the decimal point, for rounding
    point scores in \code{print.nomogram}.  Default is to round to the nearest
    whole number of points.
  }

  \item{lplabel}{
    label for linear predictor axis.  Default is \code{"Linear Predictor"}.
  }
  \item{fun.side}{
    a vector or list of vectors of \code{side} parameters for the \code{axis} function
    for labeling function values.
    Values may be 1 to position a tick mark label below the axis (the default),
    or 3 for above the axis.  If for example an axis has 5 tick mark labels
    and the second and third will run into each other, specify
    \code{fun.side=c(1,1,3,1,1)} (assuming only one function is specified as \code{fun}).
  }
  \item{col.conf}{
    colors corresponding to \code{conf.int}.
  }
  \item{conf.space}{
    a 2-element vector with the vertical range within which to draw
    confidence bars, in units of 1=spacing between main bars.  Four heights
    are used within this range (8 for the linear predictor if more than
    16 unique values were evaluated), cycling them among separate confidence
    intervals to reduce overlapping.
  }
  \item{label.every}{
    Specify \code{label.every = i} to label on every \code{i}th tick mark.
  }
  \item{force.label}{
    set to \code{TRUE} to force every tick mark intended to be labeled to have
    a label plotted (whether the labels run into each other or not)
  }
  \item{xfrac}{
    fraction of horizontal plot to set aside for axis titles
  }
  \item{cex.axis}{
    character size for tick mark labels
  }
  \item{cex.var}{
    character size for axis titles (variable names)
  }
  \item{col.grid}{
    If left unspecified, no vertical reference lines are drawn.  Specify a
    vector of length one (to use the same color for both minor and major
    reference lines) or two (corresponding to the color for the major and
    minor divisions, respectively) containing colors, to cause vertical reference
    lines to the top points scale to be drawn.  For R, a good choice is
    \code{col.grid = gray(c(0.8, 0.95))}.
  }
  \item{varname.label}{
    In constructing axis titles for interactions, the default is to add
    \code{(interacting.varname = level)} on the right.  Specify
    \code{varname.label = FALSE}
    to instead use \code{"(level)"}.
  }
  \item{varname.label.sep}{
    If \code{varname.label = TRUE}, you can change the separator to something other than
    \code{=} by specifying this parameter.
  }
  \item{ia.space}{
    When multiple axes are draw for levels of interacting factors, the
    default is to group combinations related to a main effect.  This is
    done by spacing the axes for the second to last of these 
    within a group only
    0.7 (by default) of the way down as compared with normal space of 1 unit.
  }
  \item{tck}{
    see \code{tck} under \code{\link{par}}
  }
  \item{tcl}{length of tick marks in nomogram}
  \item{lmgp}{
    spacing between numeric axis labels and axis (see \code{\link{par}}
    for \code{mgp})
  }
  \item{naxes}{
    maximum number of axes to allow on one plot.  If the nomogram requires more
    than one \dQuote{page}, the \dQuote{Points} axis will be repeated at
    the top of each page when necessary.
  }
  \item{points.label}{
    a character string giving the axis label for the points scale
  }
  \item{total.points.label}{
    a character string giving the axis label for the total points scale
  }
  \item{total.sep.page}{
    set to \code{TRUE} to force the total points and later axes to be placed on a
    separate page
  }
  \item{total.fun}{
    a user-provided function that will be executed before the total points
    axis is drawn.  Default is not to execute a function.  This is useful e.g.
    when \code{total.sep.page = TRUE} and you wish to use \code{locator} to find the
    coordinates for positioning an abbreviation legend before it's too late
    and a new page is started (i.e., \code{total.fun = function() print(locator(1))}).
  }
  \item{cap.labels}{logical: should the factor labels have their first
    letter capitalized?}

  
  \item{object}{
    the result returned from \code{nomogram}
  }
  \item{which}{
    a character string giving the name of a variable for which to draw a
    legend with abbreviations of factor levels
  }
  \item{y}{
    y-coordinate to pass to the \code{legend} function.  This is the upper left
    corner of the legend box.  You can omit \code{y} if \code{x} is a list with
    named elements \code{x} and \code{y}.  To use the mouse to locate the legend,
    specify \code{locator(1)} for \code{x}.  For \code{print}, \code{x} is
    the result of \code{nomogram}.
  }
  \item{ncol}{
    the number of columns to form in drawing the legend.
  }
}
\value{
  a list of class \code{"nomogram"} that contains information used in plotting
  the axes.  If you specified \code{abbrev = TRUE}, a list called \code{abbrev} is also
  returned that gives the abbreviations used for tick mark labels, if any.  
  This list is useful for
  making legends and is used by \code{legend.nomabbrev} (see the last example).
  The returned list also has components called \code{total.points}, \code{lp},
  and the function axis names.  These components have components
  \code{x} (\code{at} argument vector given to \code{axis}), \code{y} (\code{pos} for \code{axis}),
  and \code{x.real}, the x-coordinates appearing on tick mark labels.
  An often useful result is stored in the list of data for each axis variable,
  namely the exact number of points that correspond to each tick mark on
  that variable's axis.
}
\details{
  A variable is considered to be discrete if it is categorical or ordered
  or if \code{\link{datadist}} stored \code{values} for it (meaning it
    had \code{<11} unique values).
  A variable is said to be indirectly related to another variable if
  the two are related by some interaction.  For example, if a model
  has variables a, b, c, d, and the interactions are a:c and c:d,
  variable d is indirectly related to variable a.  The complete list
  of variables related to a is c, d.  If an axis is made for variable a,
  several axes will actually be drawn, one for each combination of c
  and d specified in \code{interact}.

  Note that with a caliper, it is easy to continually add point scores
  for individual predictors, and then to place the caliper on the upper
  \dQuote{Points} axis (with extrapolation if needed).  Then transfer these
  points to the
  \dQuote{Total Points} axis.  In this way, points can be added without
  writing them down.

  Confidence limits for an individual predictor score are really confidence
  limits for the entire linear predictor, with other predictors set to
  adjustment values.  If \code{lp = TRUE}, all confidence bars for all linear
  predictor values evaluated are drawn.  The extent to which multiple
  confidence bars of differing widths appear at the same linear predictor
  value means that precision depended on how the linear predictor was
  arrived at (e.g., a certain value may be realized from a setting of
  a certain predictor that was associated with a large standard error
  on the regression coefficients for that predictor).

  On occasion, you may want to reverse the regression coefficients of a model
  to make the \dQuote{points} scales reverse direction.  For parametric survival
  models, which are stated in terms of increasing regression effects meaning
  longer survival (the opposite of a Cox model), just do something like
  \code{fit$coefficients <- -fit$coefficients} before invoking \code{nomogram}, 
  and if you add function axes, negate the function
  arguments.  For the Cox model, you also need to negate \code{fit$center}.
  If you omit \code{lp.at}, also negate \code{fit$linear.predictors}.
}
\author{
  Frank Harrell\cr
  Department of Biostatistics\cr
  Vanderbilt University\cr
  \email{f.harrell@vanderbilt.edu}
}
\references{
  Banks J: Nomograms. Encylopedia of Statistical Sciences, Vol 6.
  Editors: S Kotz and NL Johnson.  New York: Wiley; 1985.


  Lubsen J, Pool J, van der Does, E: A practical device for the application
  of a diagnostic or prognostic function.  Meth. Inform. Med. 17:127--129;
  1978.

  Wikipedia: Nomogram, \url{http://en.wikipedia.org/wiki/Nomogram}.
}
\seealso{
  \code{\link{rms}}, \code{\link{plot.Predict}},
	\code{\link{ggplot.Predict}},	\code{\link{plot.summary.rms}},
	\code{\link{axis}}, \code{\link{pretty}}, \code{\link{approx}}, 
  \code{\link{latexrms}}, \code{\link{rmsMisc}}
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
d <- data.frame(age = rnorm(n, 50, 10),
                blood.pressure = rnorm(n, 120, 15),
                cholesterol = rnorm(n, 200, 25),
                sex = factor(sample(c('female','male'), n,TRUE)))

# Specify population model for log odds that Y=1
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
d <- upData(d,
  L = .4*(sex=='male') + .045*(age-50) +
       (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male')),
  y = ifelse(runif(n) < plogis(L), 1, 0))

ddist <- datadist(d); options(datadist='ddist')


f <- lrm(y ~ lsp(age,50) + sex * rcs(cholesterol, 4) + blood.pressure,
         data=d)
nom <- nomogram(f, fun=function(x)1/(1+exp(-x)),  # or fun=plogis
    fun.at=c(.001,.01,.05,seq(.1,.9,by=.1),.95,.99,.999),
    funlabel="Risk of Death")
#Instead of fun.at, could have specified fun.lp.at=logit of
#sequence above - faster and slightly more accurate
plot(nom, xfrac=.45)
print(nom)
nom <- nomogram(f, age=seq(10,90,by=10))
plot(nom, xfrac=.45)
g <- lrm(y ~ sex + rcs(age, 3) * rcs(cholesterol, 3), data=d)
nom <- nomogram(g, interact=list(age=c(20,40,60)), 
                conf.int=c(.7,.9,.95))
plot(nom, col.conf=c(1,.5,.2), naxes=7)

w <- upData(d,
            cens = 15 * runif(n),
            h = .02 * exp(.04 * (age - 50) + .8 * (sex == 'Female')),
            d.time = -log(runif(n)) / h,
            death = ifelse(d.time <= cens, 1, 0),
            d.time = pmin(d.time, cens))


f <- psm(Surv(d.time,death) ~ sex * age, data=w, dist='lognormal')
med  <- Quantile(f)
surv <- Survival(f)  # This would also work if f was from cph
plot(nomogram(f, fun=function(x) med(lp=x), funlabel="Median Survival Time"))
nom <- nomogram(f, fun=list(function(x) surv(3, x),
                            function(x) surv(6, x)),
            funlabel=c("3-Month Survival Probability", 
                       "6-month Survival Probability"))
plot(nom, xfrac=.7)

\dontrun{
nom <- nomogram(fit.with.categorical.predictors, abbrev=TRUE, minlength=1)
nom$x1$points   # print points assigned to each level of x1 for its axis
#Add legend for abbreviations for category levels
abb <- attr(nom, 'info')$abbrev$treatment
legend(locator(1), abb$full, pch=paste(abb$abbrev,collapse=''), 
       ncol=2, bty='n')  # this only works for 1-letter abbreviations
#Or use the legend.nomabbrev function:
legend.nomabbrev(nom, 'treatment', locator(1), ncol=2, bty='n')
}


#Make a nomogram with axes predicting probabilities Y>=j for all j=1-3
#in an ordinal logistic model, where Y=0,1,2,3
w <- upData(w, Y = ifelse(y==0, 0, sample(1:3, length(y), TRUE)))
g <- lrm(Y ~ age+rcs(cholesterol,4) * sex, data=w)
fun2 <- function(x) plogis(x-g$coef[1]+g$coef[2])
fun3 <- function(x) plogis(x-g$coef[1]+g$coef[3])
f <- Newlabels(g, c(age='Age in Years'))  
#see Design.Misc, which also has Newlevels to change 
#labels for levels of categorical variables
g <- nomogram(f, fun=list('Prob Y>=1'=plogis, 'Prob Y>=2'=fun2, 
                     'Prob Y=3'=fun3), 
         fun.at=c(.01,.05,seq(.1,.9,by=.1),.95,.99))
plot(g, lmgp=.2, cex.axis=.6)
options(datadist=NULL)
}
\keyword{models}
\keyword{regression}
\keyword{hplot}


