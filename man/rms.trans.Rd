\name{rms.trans}
\alias{rms.trans}
\alias{asis}
\alias{pol}
\alias{lsp}
\alias{rcs}
\alias{catg}
\alias{scored}
\alias{strat}
\alias{matrx}
\alias{gTrans}
\alias{\%ia\%}
\alias{makepredictcall.rms}
\title{rms Special Transformation Functions}
\description{
This is a series of functions (\code{asis}, \code{pol}, \code{lsp},
\code{rcs}, \code{catg}, \code{scored}, \code{strat}, \code{matrx},
\code{gTrans}, and
\code{\%ia\%}) that set up special attributes  (such as
knots and nonlinear term indicators) that are carried through to fits
(using for example \code{lrm},\code{cph}, \code{ols},
\code{psm}). \code{anova.rms}, \code{summary.rms}, \code{Predict},
\code{survplot}, \code{fastbw}, \code{validate}, \code{specs},
\code{which.influence}, \code{nomogram} and \code{latex.rms} use these 
attributes to automate certain analyses (e.g., automatic tests of linearity
for each predictor are done by \code{anova.rms}).  Many of the functions
are called implicitly.  Some S functions such as \code{ns} derive data-dependent
transformations that are not always "remembered" when predicted values are
later computed, so the predictions may be incorrect. The functions listed
here solve that problem when used in the \code{rms} context.

\code{asis} is the identity transformation, \code{pol} is an ordinary
(non-orthogonal) polynomial, \code{rcs} is a linear tail-restricted
cubic spline function (natural spline, for which the
\code{rcspline.eval} function generates the design matrix, the
presence of system option \code{rcspc} causes \code{rcspline.eval} to be
invoked with \code{pc=TRUE}, and the presence of system option \code{fractied}
causes this value to be passed to \code{rcspline.eval} as the \code{fractied}
argument), \code{catg} is for a categorical variable,
\code{scored} is for an ordered categorical variable, \code{strat} is
for a stratification factor in a Cox model, \code{matrx} is for a matrix
predictor, and \code{\%ia\%} represents restricted interactions in which
products involving nonlinear effects on both variables are not included
in the model.  \code{asis, catg, scored, matrx} are seldom invoked
explicitly by the user (only to specify \code{label} or \code{name},
usually).

\code{gTrans} is a general multiple-parameter transformation function.
It can be used to specify new polynomial bases, smooth relationships
with a discontinuity at one or more values of \code{x}, grouped
categorical variables, e.g., a categorical variable with 5 levels where
you want to combine two of the levels to spend only 3 degrees of freedom in
all but see plots of predicted values where the two combined categories
are kept separate but will have equal effect estimates.  The first
argument to \code{gTrans} is a regular numeric, character, or factor
variable.  The next argument is a function that transforms a vector into
a matrix.  If the basis functions are to include a linear term it is up
too the user to include the original \code{x} as one of the columns.
Column names are assigned automaticall, but any column names specified
by the user will override the default name.  If you want to signal which
terms correspond to linear and which correspond to nonlinear effects for
the purpose of running \code{anova.rms}, add an integer vector attribute
\code{nonlinear} to the resulting matrix.  This vector specifies the
column numbers corresponding to nonlinear effects.  The default is to assume a column
is a linear effect.  The \code{parms} attribute stored with a
\code{gTrans} result a character vector version of the function, so as
to not waste space carrying along any environment information.  If you
will be using the \code{latex} method for typesetting the fitted model,
you must include a \code{tex} attribute also in the produced matrix.
This must be a function of a single character string argument (that will
ultimately contain the name of the predictor in LaTeX notation) and must
produce a vector of LaTeX character strings.  See
\url{https://hbiostat.org/R/rms/gtrans.html} for several examples of the
use of \code{gTrans} including the use of \code{nonlinear} and
\code{tex}.

A \code{makepredictcall} method is defined so that usage of the
transformation functions outside of \code{rms} fitting functions will
work for getting predicted values.  Thanks to Therry Therneau for the code.

In the list below, functions \code{asis} through \code{gTrans} can have
arguments \code{x, parms, label, name} except that \code{parms} does not
apply to \code{asis, matrx, strat}.
}
\usage{
asis(\dots)
matrx(\dots)
pol(\dots)
lsp(\dots)
rcs(\dots)
catg(\dots)
scored(\dots)
strat(\dots)
gTrans(\dots)
x1 \%ia\% x2
\method{makepredictcall}{rms}(var, call)
}
\arguments{
  \item{\dots}{
	The arguments \dots above contain the following.
	\describe{
	  \item{\code{x}}{a predictor variable (or a function of one).  If
	  you specify e.g. \code{pol(pmin(age,10),3)}, a cubic polynomial
	  will be fitted in \code{pmin(age,10)} (\code{pmin} is the S vector
	  element--by--element function). The predictor will be labeled
	  \code{age} in the output, and plots with have \code{age} in its
	  original units on the axes. If you use a function such as
	  \code{pmin}, the predictor is taken as the first argument, and
	  other arguments must be defined in the frame in effect when
	  predicted values, etc., are computed.}
	  \item{\code{parms}}{parameters of transformation (e.g. number or
	  location of knots). For \code{pol} the argument is the order of
	  the polynomial, e.g. \code{2} for quadratic (the usual
	  default). For \code{lsp} it is a vector of knot locations
	  (\code{lsp} will not estimate knot locations).  For \code{rcs} it
	  is the number of knots (if scalar), or vector of knot locations
	  (if \code{>2} elements).  The default number is the \code{nknots}
	  system option if \code{parms} is not given.  If the number of
	  knots is given, locations are computed for that number of knots.
	  If system option \code{rcspc} is \code{TRUE} the \code{parms}
	  vector has an attribute defining the principal components
	  transformation parameters.  For \code{catg}, \code{parms} is the
	  category labels (not needed if variable is an S category or factor
	  variable). If omitted, \code{catg} will use \code{unique(x)}, or
	  \code{levels(x)} if \code{x} is a \code{category} or a
	  \code{factor}.  For \code{scored}, \code{parms} is a vector of
	  unique values of variable (uses \code{unique(x)} by default).
	  This is not needed if \code{x} is an S \code{ordered} variable.
	  For \code{strat}, \code{parms} is the category labels (not needed
	  if variable is an S category variable). If omitted, will use
	  \code{unique(x)}, or \code{levels(x)} if \code{x} is
	  \code{category} or \code{factor}. \code{parms} is not used for
	  \code{matrix}.}
	  \item{\code{label}}{label of predictor for plotting (default =
	  \code{"label"} attribute or variable name)}
	  \item{\code{name}}{Name to use for predictor in model. Default is
	  name of argument to function.}
	}
}
\item{x1,x2}{two continuous variables for which to form a
  non-doubly-nonlinear interaction}
\item{var}{a model term passed from a (usually non-\code{rms}) function}
\item{call}{call object for a model term}
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link[Hmisc]{rcspline.eval}},
\code{\link[Hmisc]{rcspline.restate}}, \code{\link{rms}},
\code{\link{cph}}, \code{\link{lrm}}, \code{\link{ols}},
\code{\link{datadist}}, \code{\link[stats]{makepredictcall}}
}
\examples{
\dontrun{
options(knots=4, poly.degree=2)
# To get the old behavior of rcspline.eval knot placement (which didnt' handle
# clumping at the lowest or highest value of the predictor very well):
# options(fractied = 1.0)   # see rcspline.eval for details
country <- factor(country.codes)
blood.pressure <- cbind(sbp=systolic.bp, dbp=diastolic.bp)
fit <- lrm(Y ~ sqrt(x1)*rcs(x2) + rcs(x3,c(5,10,15)) + 
       lsp(x4,c(10,20)) + country + blood.pressure + poly(age,2))
# sqrt(x1) is an implicit asis variable, but limits of x1, not sqrt(x1)
#       are used for later plotting and effect estimation
# x2 fitted with restricted cubic spline with 4 default knots
# x3 fitted with r.c.s. with 3 specified knots
# x4 fitted with linear spline with 2 specified knots
# country is an implied catg variable
# blood.pressure is an implied matrx variable
# since poly is not an rms function (pol is), it creates a
#       matrx type variable with no automatic linearity testing
#       or plotting
f1 <- lrm(y ~ rcs(x1) + rcs(x2) + rcs(x1) \%ia\% rcs(x2))
# \%ia\% restricts interactions. Here it removes terms nonlinear in
# both x1 and x2
f2 <- lrm(y ~ rcs(x1) + rcs(x2) + x1 \%ia\% rcs(x2))
# interaction linear in x1
f3 <- lrm(y ~ rcs(x1) + rcs(x2) + x1 \%ia\% x2)
# simple product interaction (doubly linear)
# Use x1 \%ia\% x2 instead of x1:x2 because x1 \%ia\% x2 triggers
# anova to pool x1*x2 term into x1 terms to test total effect
# of x1
#
# Examples of gTrans
#
# Linear relationship with a discontinuity at zero:
ldisc <- function(x) {z <- cbind(x == 0, x); attr(z, 'nonlinear') <- 1; z}
gTrans(x, ldisc)
# Duplicate pol(x, 2):
pol2 <- function(x) {z <- cbind(x, x^2); attr(z, 'nonlinear') <- 2; z}
gTrans(x, pol2)
# Linear spline with a knot at x=10 with the new slope taking effect
# until x=20 and the spline turning flat at that point but with a
# discontinuous vertical shift
# tex is only needed if you will be using latex(fit)
dspl <- function(x) {
  z <- cbind(x, pmax(pmin(x, 20) - 10, 0), x > 20)
  attr(z, 'nonlinear') <- 2:3
  attr(z, 'tex') <- function(x) sprintf(c('%s', '(\\min(%s, 20) - 10)_{+}',
                                          '[%s > 20]'), x)
  z }
gTrans(x, dspl)
}
}
\keyword{models}
\keyword{regression}
\keyword{math}
\keyword{manip}
\keyword{methods}
\keyword{survival}
\keyword{smooth}
\concept{logistic regression model}
\concept{transformation}
