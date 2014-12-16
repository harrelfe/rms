\name{Gls}
\alias{Gls}
\alias{print.Gls}
\title{Fit Linear Model Using Generalized Least Squares}
\usage{
Gls(model, data, correlation, weights, subset, method, na.action=na.omit,
    control, verbose, B=0, dupCluster=FALSE, pr=FALSE, x=FALSE)

\method{print}{Gls}(x, digits=4, coefs=TRUE, latex=FALSE, title, \dots)
}
\arguments{
 \item{model}{a two-sided linear formula object describing the
    model, with the response on the left of a \code{~} operator and the
    terms, separated by \code{+} operators, on the right.}
 \item{data}{an optional data frame containing the variables named in
   \code{model}, \code{correlation}, \code{weights}, and
   \code{subset}. By default the variables are taken from the
   environment from which \code{gls} is called.}
 \item{correlation}{an optional \code{corStruct} object describing the
   within-group correlation structure. See the documentation of
   \code{corClasses} for a description of the available \code{corStruct}
   classes. If a grouping variable is to be used, it must be specified in
   the \code{form} argument to the \code{corStruct}
   constructor. Defaults to \code{NULL}, corresponding to uncorrelated 
   errors.}  
 \item{weights}{an optional \code{varFunc} object or one-sided formula
   describing the within-group heteroscedasticity structure. If given as
   a formula, it is used as the argument to \code{varFixed},
   corresponding to fixed variance weights. See the documentation on
   \code{varClasses} for a description of the available \code{varFunc}
   classes. Defaults to \code{NULL}, corresponding to homoscesdatic
   errors.} 
 \item{subset}{an optional expression indicating which subset of the rows of
   \code{data} should  be  used in the fit. This can be a logical
   vector, or a numeric vector indicating which observation numbers are
   to be included, or a  character  vector of the row names to be
   included.  All observations are included by default.}
 \item{method}{a character string.  If \code{"REML"} the model is fit by
   maximizing the restricted log-likelihood.  If \code{"ML"} the
   log-likelihood is maximized.  Defaults to \code{"REML"}.}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.omit}) results
   in deletion of observations having any of the variables of interest missing.}
 \item{control}{a list of control values for the estimation algorithm to
   replace the default values returned by the function \code{glsControl}.
   Defaults to an empty list.}
 \item{verbose}{an optional logical value. If \code{TRUE} information on
   the evolution of the iterative algorithm is printed. Default is
   \code{FALSE}.}
 \item{B}{number of bootstrap resamples to fit and store, default is
   none}
 \item{dupCluster}{set to \code{TRUE} to have \code{Gls} when
   bootstrapping to consider multiply-sampled clusters as if they were
   one large cluster when fitting using the \code{gls} algorithm}
 \item{pr}{set to \code{TRUE} to show progress of bootstrap resampling}
 \item{x}{for \code{Gls} set to \code{TRUE} to store the design matrix
   in the fit object; otherwise the result of \code{Gls}}
 \item{digits}{number of significant digits to print}
 \item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
  of model coefficients, standard errors, etc.  Specify \code{coefs=n}
  to print only the first \code{n} regression coefficients in the
  model.}
\item{latex}{a logical value indicating whether information should be
  formatted as plain text or as LaTeX markup}
\item{title}{a character string title to be passed to \code{prModFit}}
 \item{\dots}{ignored}
}
\description{
  This function fits a linear model using generalized least
  squares. The errors are allowed to be correlated and/or have unequal
  variances.  \code{Gls} is a slightly enhanced version of the
  Pinheiro and Bates \code{Gls} function in the \code{nlme} package to
  make it easy to use with the rms package and to implement cluster
  bootstrapping (primarily for nonparametric estimates of the
  variance-covariance matrix of the parameter estimates and for
  nonparametric confidence limits of correlation parameters).
}
\value{
  an object of classes \code{Gls}, \code{rms}, and \code{gls}
  representing the linear model
  fit. Generic functions such as \code{print}, \code{plot},
	\code{ggplot}, and \code{summary} have methods to show the results of
	the fit. See 
  \code{glsObject} for the components of the fit. The functions
  \code{resid}, \code{coef}, and \code{fitted} can be used to extract
  some of its components.  \code{Gls} returns the following components
  not returned by \code{gls}: \code{Design}, \code{assign},
  \code{formula} (see arguments), \code{B} (see
  arguments), \code{bootCoef} (matrix of \code{B} bootstrapped
  coefficients), \code{boot.Corr} (vector of bootstrapped correlation
  parameters), \code{Nboot} (vector of total sample size used in each
  bootstrap (may vary if have unbalanced clusters), and \code{var}
  (sample variance-covariance matrix of bootstrapped coefficients).  The
  \eqn{g}-index is also stored in the returned object under the name
  \code{"g"}.
}
\references{
  Pinheiro J, Bates D (2000): Mixed effects models in S and S-Plus.  New
  York: Springer-Verlag.
}
\details{
  The \code{\link[Hmisc]{na.delete}} function will not work with
  \code{Gls} due to some nuance in the \code{model.frame.default}
  function.  This probably relates to \code{na.delete} storing extra
  information in the \code{"na.action"} attribute of the returned data
  frame.
}
\author{Jose Pinheiro,
  Douglas Bates \email{bates@stat.wisc.edu},
	Saikat DebRoy,
	Deepayan Sarkar,
	R-core \email{R-core@R-project.org},
  Frank Harrell \email{f.harrell@vanderbilt.edu},
  Patrick Aboyoun
} 
\seealso{
  \code{\link[nlme]{gls}}
  \code{\link[nlme]{glsControl}}, \code{\link[nlme]{glsObject}},
  \code{\link{varFunc}}, \code{\link[nlme]{corClasses}},
  \code{\link[nlme]{varClasses}}, \code{\link{GiniMd}},
	\code{\link{prModFit}}, \code{\link{logLik.Gls}}
}
\examples{
\dontrun{
ns  <- 20  # no. subjects
nt  <- 10  # no. time points/subject
B   <- 10  # no. bootstrap resamples
           # usually do 100 for variances, 1000 for nonparametric CLs
rho <- .5  # AR(1) correlation parameter
V <- matrix(0, nrow=nt, ncol=nt)
V <- rho^abs(row(V)-col(V))   # per-subject correlation/covariance matrix

d <- expand.grid(tim=1:nt, id=1:ns)
d$trt <- factor(ifelse(d$id <= ns/2, 'a', 'b'))
true.beta <- c(Intercept=0,tim=.1,'tim^2'=0,'trt=b'=1)
d$ey  <- true.beta['Intercept'] + true.beta['tim']*d$tim +
  true.beta['tim^2']*(d$tim^2) +  true.beta['trt=b']*(d$trt=='b')
set.seed(13)
library(MASS)   # needed for mvrnorm
d$y <- d$ey + as.vector(t(mvrnorm(n=ns, mu=rep(0,nt), Sigma=V)))

dd <- datadist(d); options(datadist='dd')
f <- Gls(y ~ pol(tim,2) + trt, correlation=corCAR1(form= ~tim | id),
         data=d, B=B)
f
AIC(f)
f$var      # bootstrap variances
f$varBeta  # original variances
summary(f)
anova(f)
ggplot(Predict(f, tim, trt))
# v <- Variogram(f, form=~tim|id, data=d)
nlme:::summary.gls(f)$tTable   # print matrix of estimates etc.

options(datadist=NULL)
}
}
\keyword{models}
