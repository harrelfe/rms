\name{rmsMisc}
\alias{rmsMisc}
\alias{calibrate.rms}
\alias{vcov.rms}
\alias{vcov.cph}
\alias{vcov.Glm}
\alias{vcov.Gls}
\alias{vcov.lrm}
\alias{vcov.ols}
\alias{vcov.orm}
\alias{vcov.psm}
\alias{oos.loglik}
\alias{oos.loglik.ols}
\alias{oos.loglik.lrm}
\alias{oos.loglik.cph}
\alias{oos.loglik.psm}
\alias{oos.loglik.Glm}
\alias{Getlim}
\alias{Getlimi}
\alias{related.predictors}
\alias{interactions.containing}
\alias{combineRelatedPredictors}
\alias{param.order}
\alias{Penalty.matrix}
\alias{Penalty.setup}
\alias{logLik.Gls}
\alias{logLik.ols}
\alias{logLik.rms}
\alias{AIC.rms}
\alias{nobs.rms}
\alias{lrtest}
\alias{univarLR}
\alias{Newlabels}
\alias{Newlevels}
\alias{Newlabels.rms}
\alias{Newlevels.rms}
\alias{rmsArgs}
\alias{print.rms}
\alias{print.lrtest}
\alias{survest.rms}
\alias{prModFit}
\alias{prStats}
\alias{reVector}
\alias{formatNP}
\alias{latex.naprint.delete}
\alias{removeFormulaTerms}
\title{Miscellaneous Design Attributes and Utility Functions}
\description{
These functions are used internally to \code{anova.rms},
\code{fastbw}, etc., to retrieve various attributes of a design.  These
functions allow some fitting functions not in the \code{rms} series
(e.g,, \code{lm}, \code{glm}) to be used with \code{rms.Design},
\code{fastbw}, and similar functions. 


For \code{vcov}, there are several functions.  The method for \code{orm}
fits is a bit different because the covariance matrix stored in the fit
object only deals with the middle intercept.  See the \code{intercepts}
argument for more options.  There is a method for \code{lrm} that also
allows non-default intercept(s) to be selected (default is first).

The \code{oos.loglik} function for
each type of model implemented computes the -2 log likelihood for
out-of-sample data (i.e., data not necessarily used to fit the model)
evaluated at the parameter estimates from a model fit.  Vectors for the
model's linear predictors and response variable must be given.
\code{oos.loglik} is used primarily by \code{bootcov}.


The \code{Getlim} function retrieves distribution summaries
from the fit or from a \code{datadist} object.  It handles getting summaries
from both sources to fill in characteristics for variables that were not
defined during the model fit.  \code{Getlimi} returns the summary
for an individual model variable.  


The \code{related.predictors} function
returns a list containing variable numbers that are directly or
indirectly related to each predictor.  The \code{interactions.containing}
function returns indexes of interaction effects containing a given
predictor.  The \code{param.order} function returns a vector of logical
indicators for whether parameters are associated with certain types of
effects (nonlinear, interaction, nonlinear interaction).
\code{combineRelatedPredictors} creates of list of inter-connected main
effects and interations for use with \code{predictrms} with
\code{type='ccterms'} (useful for \code{gIndex}).


The \code{Penalty.matrix} function builds a default penalty matrix for
non-intercept term(s) for use in penalized maximum likelihood
estimation.  The \code{Penalty.setup} function takes a constant or list
describing penalty factors for each type of term in the model and
generates the proper vector of penalty multipliers for the current model.

\code{logLik.rms} returns the maximized log likelihood for the model,
whereas \code{AIC.rms} returns the AIC.  The latter function has an
optional argument for computing AIC on a "chi-square" scale (model
likelihood ratio chi-square minus twice the regression degrees of
freedom.  \code{logLik.ols} handles the case for \code{ols}, just by
invoking \code{logLik.lm} in the \code{stats} package.
\code{logLik.Gls} is also defined.

\code{nobs.rms} returns the number of observations used in the fit.

The \code{lrtest} function does likelihood ratio tests for
two nested models, from fits that have \code{stats} components with
\code{"Model L.R."} values.  For models such as \code{psm, survreg, ols, lm} which have
scale parameters, it is assumed that scale parameter for the smaller model
is fixed at the estimate from the larger model (see the example).

\code{univarLR} takes a multivariable model fit object from
\code{rms} and re-fits a sequence of models containing one predictor
at a time.  It prints a table of likelihood ratio \eqn{chi^2} statistics
from these fits.

The \code{Newlabels} function is used to override the variable labels in a
fit object.  Likewise, \code{Newlevels} can be used to create a new fit object
with levels of categorical predictors changed.  These two functions are
 especially useful when constructing nomograms.

\code{rmsArgs} handles \dots arguments to functions such as
\code{Predict}, \code{summary.rms}, \code{nomogram} so that variables to
vary may be specified without values (after an equals sign).

\code{prModFit} is the workhorse for the \code{print} methods for
highest-level \code{rms} model fitting functions, handling both regular
and LaTeX printing, the latter resulting in LaTeX code written to the
terminal, automatically ready for \code{Sweave}.  The work of printing
summary statistics is done by \code{prStats}, which uses the Hmisc
\code{print.char.matrix} function to print overall model statistics if
\code{latex=FALSE}, otherwise it generates customized LaTeX code.  The
LaTeX longtable and epic packages must be in effect to use these
LaTeX functions.

\code{reVector} allows one to rename a subset of a named vector,
ignoring the previous names and not concatenating them as \R does.  It
also removes (by default) elements that are \code{NA}, as when an
optional named element is fetched that doesn't exist.

\code{formatNP} is a function to format a vector of numerics.  If
\code{digits} is specified, \code{formatNP} will make sure that the
formatted representation has \code{digits} positions to the right of the
decimal place.  If \code{latex=TRUE} it will translate any scientific
notation to LaTeX math form.  If \code{pvalue=TRUE}, it will replace
formatted values with "< 0.0001" (if \code{digits=4}).

\code{latex.naprint.delete} will, if appropriate, use LaTeX to draw a
dot chart of frequency of variable \code{NA}s related to model fits.

\code{removeFormulaTerms} removes one or more terms from a model
formula, using strictly character manipulation.  This handles problems
such as \code{[.terms} removing \code{offset()} if you subset on
anything.  The function can also be used to remove the dependent
variable(s) from the formula.
}
\usage{
\method{vcov}{rms}(object, regcoef.only=TRUE, intercepts='all', \dots)
\method{vcov}{cph}(object, regcoef.only=TRUE, \dots)
\method{vcov}{Glm}(object, regcoef.only=TRUE, intercepts='all', \dots)
\method{vcov}{Gls}(object, intercepts='all', \dots)
\method{vcov}{lrm}(object, regcoef.only=TRUE, intercepts='all', \dots)
\method{vcov}{ols}(object, regcoef.only=TRUE, \dots)
\method{vcov}{orm}(object, regcoef.only=TRUE, intercepts='mid', \dots)
\method{vcov}{psm}(object, regcoef.only=TRUE, \dots)

oos.loglik(fit, \dots)

\method{oos.loglik}{ols}(fit, lp, y, \dots)
\method{oos.loglik}{lrm}(fit, lp, y, \dots)
\method{oos.loglik}{cph}(fit, lp, y, \dots)
\method{oos.loglik}{psm}(fit, lp, y, \dots)
\method{oos.loglik}{Glm}(fit, lp, y, \dots)

Getlim(at, allow.null=FALSE, need.all=TRUE)
Getlimi(name, Limval, need.all=TRUE)

related.predictors(at, type=c("all","direct"))
interactions.containing(at, pred)
combineRelatedPredictors(at)
param.order(at, term.order)

Penalty.matrix(at, X)
Penalty.setup(at, penalty)

\method{logLik}{Gls}(object, \dots)
\method{logLik}{ols}(object, \dots)
\method{logLik}{rms}(object, \dots)
\method{AIC}{rms}(object, \dots, k=2, type=c('loglik', 'chisq'))
\method{nobs}{rms}(object, \dots)

lrtest(fit1, fit2)
\method{print}{lrtest}(x, \dots)

univarLR(fit)

Newlabels(fit, \dots)
Newlevels(fit, \dots)
\method{Newlabels}{rms}(fit, labels, \dots)
\method{Newlevels}{rms}(fit, levels, \dots)

prModFit(x, title, w, digits=4, coefs=TRUE, latex=FALSE, rmarkdown=FALSE,
         lines.page=40, long=TRUE, needspace, \dots)

prStats(labels, w, latex=FALSE, file="", append=TRUE)

reVector(\dots, na.rm=TRUE)

formatNP(x, digits=NULL, pvalue=FALSE, latex=FALSE)

\method{latex}{naprint.delete}(object, file="", append=TRUE, \dots)

removeFormulaTerms(form, which=NULL, delete.response=FALSE)
}
\arguments{
\item{fit}{result of a fitting function}
\item{object}{result of a fitting function}
\item{regcoef.only}{For fits such as parametric survival models
  which have a final row and column of the covariance matrix for a
  non-regression parameter such as a log(scale) parameter, setting
  \code{regcoef.only=TRUE} causes only the first 
  \code{p} rows and columns of the covariance matrix to be returned,
  where \code{p} is the length of \code{object$coef}.
}
\item{intercepts}{set to \code{"none"} to omit any rows and columns
	related to intercepts.  Set to an integer scalar
	or vector to include particular intercept elements.  Set to
	\code{'all'} to include all intercepts, or for \code{orm} to
	\code{"mid"} to use the default for \code{orm}.  The default is to use the
	first for \code{lrm} and the median intercept for \code{orm}.
	}
\item{at}{
\code{Design} element of a fit
}
\item{pred}{
index of a predictor variable (main effect)
}
\item{fit1}{
}
\item{fit2}{
fit objects from \code{lrm,ols,psm,cph} etc.  It doesn't matter which
fit object is the sub-model.
}
\item{lp}{
linear predictor vector for \code{oos.loglik}.  For proportional odds
ordinal logistic models, this should have used the first intercept
only.  If \code{lp} and \code{y} are omitted, the -2 log likelihood for the
original fit are returned.
}
\item{y}{
values of a new vector of responses passed to \code{oos.loglik}.
}
\item{name}{
the name of a variable in the model
}
\item{Limval}{
an object returned by \code{Getlim}
}
\item{allow.null}{
prevents \code{Getlim} from issuing an error message if no limits are found
in the fit or in the object pointed to by \code{options(datadist=)}
}
\item{need.all}{
set to \code{FALSE} to prevent \code{Getlim} or \code{Getlimi} from issuing an error message
if data for a variable are not found
}
\item{type}{
For \code{related.predictors}, set to \code{"direct"} to return lists of
indexes of directly related factors only (those in interactions with the
predictor).  For \code{AIC.rms}, \code{type} specifies the basis on
which to return AIC.  The default is minus twice the maximized log
likelihood plus \code{k} times the degrees of freedom counting
intercept(s).  Specify \code{type='chisq'} to get a penalized model
likelihood ratio chi-square instead.  
}
\item{term.order}{
1 for all parameters, 2 for all parameters associated with either nonlinear
or interaction effects, 3 for nonlinear effects (main or interaction),
4 for interaction effects, 5 for nonlinear interaction effects.
}
\item{X}{
a design matrix, not including columns for intercepts
}
\item{penalty}{
a vector or list specifying penalty multipliers for types of model terms
}
\item{k}{the multiplier of the degrees of freedom to be used in
  computing AIC.  The default is 2.}
\item{x}{a result of \code{lrtest}, or the result of a high-level model
  fitting function (for \code{prModFit}}
\item{labels}{
a character vector specifying new labels for variables in a fit.
To give new labels for all variables, you can specify \code{labels} of the
form \code{labels=c("Age in Years","Cholesterol")}, where the list of new labels is
assumed to be the length of all main effect-type variables in the fit and
in their original order in the model formula.  You may specify a named
vector to give new labels in random order or for a subset of the 
variables, e.g., \code{labels=c(age="Age in Years",chol="Cholesterol")}.
For \code{prStats}, is a list with major column headings, which can
themselves be vectors that are then stacked vertically.
}
\item{levels}{
a list of named vectors specifying new level labels for categorical
predictors.  This will override \code{parms} as well as \code{datadist} information
(if available) that were stored with the fit.  
}
\item{title}{a single character string used to specify an overall title
  for the regression fit, which is printed first by \code{prModFit}.
  Set to \code{""} to suppress the title}
\item{w}{For \code{prModFit}, a special list of lists, which each list
  element specifying information about a block of information to include
  in the \code{print.} output for a fit.  For \code{prStats}, \code{w}
  is a list of statistics to print, elements of which can be vectors
  that are stacked vertically.  Unnamed elements specify number of
  digits to the right of the decimal place to which to round (\code{NA}
  means use \code{format} without rounding, as with integers and
  floating point values).  Negative values of \code{digits} indicate
  that the value is a P-value to be formatted with \code{formatNP}.
  Digits are recycled as needed. 
}
\item{digits}{number of digits to the right of the decimal point, for
  formatting numeric values in printed output}
\item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
  of model coefficients, standard errors, etc.  Specify \code{coefs=n}
  to print only the first \code{n} regression coefficients in the
  model.}
\item{latex}{a logical value indicating whether information should be
  formatted as plain text or as LaTeX markup}
\item{file}{name of file to which to write model output from
	\code{print()} using \code{prStats}.  Default is the console.}
\item{append}{specify \code{append=FALSE} when using \code{file} and you
	want to start over instead of adding to an existing file.}
\item{rmarkdown}{set to \code{TRUE} to force \code{latex=TRUE} and to
	convert LaTeX code to html using Hmisc \code{html.latex} for use with
	RMarkdown, knitr, and RStudio}
\item{lines.page}{see \code{\link[Hmisc]{latex}}}
\item{long}{set to \code{FALSE} to suppress printing of formula and
  certain other model output}
\item{needspace}{optional character string to insert inside a LaTeX
  needspace macro call before the statistics table and before the
  coefficient matrix, to avoid bad page splits.  This assumes the LaTeX
  needspace style is available.  Example:
  \code{needspace='6\\baselineskip'} or \code{needspace='1.5in'}.}
\item{na.rm}{set to \code{FALSE} to keep \code{NA}s in the vector
  created by \code{reVector}}
\item{pvalue}{set to \code{TRUE} if you want values below 10 to the
  minus \code{digits} to be formatted to be less than that value}
\item{form}{a formula object}
\item{which}{a vector of one or more character strings specifying the
	names of functions that are called from a formula, e.g.,
	\code{"cluster"}.  By default no right-hand-side terms are removed.}
\item{delete.response}{set to \code{TRUE} to remove the dependent
	variable(s) from the formula}
\item{\dots}{other arguments.  For \code{reVector} this contains the
  elements being extracted.  For \code{prModFit} this information is
  passed to the \code{Hmisc latexTabular} function when a block of
  output is a vector to be formatted in LaTeX.}
}
\value{
\code{vcov} returns a variance-covariance matrix
\code{oos.loglik} returns a scalar -2 log likelihood value.
\code{Getlim} returns a list with components \code{limits} and \code{values}, either
stored in \code{fit} or retrieved from the object created by \code{datadist} and
pointed to in \code{options(datadist=)}.
\code{related.predictors} and \code{combineRelatedPredictors} return a
list of vectors, and \code{interactions.containing} 
returns a vector.  \code{param.order} returns a logical vector corresponding
to non-strata terms in the model.
\code{Penalty.matrix} returns a symmetric matrix with dimension equal to the
number of slopes in the model.  For all but categorical predictor main
effect elements, the matrix is diagonal with values equal to the variances
of the columns of \code{X}.  For segments corresponding to \code{c-1} dummy variables
for \code{c}-category predictors,  puts a \code{c-1} x \code{c-1} sub-matrix in
\code{Penalty.matrix} that is constructed so that a quadratic form with 
\code{Penalty.matrix} in the middle computes the sum of squared differences
in parameter values about the mean, including a portion for the reference
cell in which the parameter is by definition zero.
\code{Newlabels} returns a new fit object with the labels adjusted.

\code{reVector} returns a vector of named (by its arguments) elements.
\code{formatNP} returns a character vector.

\code{removeFormulaTerms} returns a formula object.
}
\seealso{
\code{\link{rms}}, \code{\link{fastbw}}, \code{\link{anova.rms}},
\code{\link{summary.lm}}, \code{\link{summary.glm}},
\code{\link{datadist}}, \code{\link{vif}}, \code{\link{bootcov}},
\code{\link[Hmisc]{latex}}, \code{\link[Hmisc]{latexTabular}},
\code{\link[Hmisc]{latexSN}}, \code{\link[Hmisc]{print.char.matrix}}
}
\examples{
\dontrun{
f <- psm(S ~ x1 + x2 + sex + race, dist='gau')
g <- psm(S ~ x1 + sex + race, dist='gau', 
         fixed=list(scale=exp(f$parms)))
lrtest(f, g)


g <- Newlabels(f, c(x2='Label for x2'))
g <- Newlevels(g, list(sex=c('Male','Female'),race=c('B','W')))
nomogram(g)
}
}
\keyword{models}
\keyword{methods}






