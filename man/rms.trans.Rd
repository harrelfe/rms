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
\alias{\%ia\%}
\title{
rms Special Transformation Functions
}
\description{
This is a series of functions (\code{asis}, \code{pol}, \code{lsp},
\code{rcs}, \code{catg}, \code{scored}, \code{strat}, \code{matrx}, and
\code{\%ia\%}) that set up special attributes  (such as
knots and nonlinear term indicators) that are carried through to fits
(using for example \code{lrm},\code{cph}, \code{ols},
\code{psm}). \code{anova.rms}, \code{summary.rms}, \code{Predict},
\code{survplot}, \code{fastbw}, \code{validate}, \code{specs},
\code{which.influence}, \code{nomogram} and \code{latex.rms} use these 
attributes to automate certain analyses (e.g., automatic tests of linearity
for each predictor are done by \code{anova.rms}). Many of the functions
are called implicitly.  Some S functions such as \code{ns} derive data-dependent
transformations that are not "remembered" when predicted values are
later computed, so the predictions will be incorrect. The functions listed
here solve that problem. 

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

In the list below, functions \code{asis} through \code{strat} can have
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
x1 \%ia\% x2
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
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link[Hmisc]{rcspline.eval}}, \code{\link[Hmisc]{rcspline.restate}}, \code{\link{rms}}, \code{\link{cph}}, \code{\link{lrm}}, \code{\link{ols}}, \code{\link{datadist}}
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
