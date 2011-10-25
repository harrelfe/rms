\name{lrm}
\alias{lrm}
\alias{print.lrm}
\title{Logistic Regression Model}

\description{
Fit binary and proportional odds ordinal
logistic regression models using maximum likelihood estimation or
penalized maximum likelihood estimation.  See \code{cr.setup} for how to
fit forward continuation ratio models with \code{lrm}.
}
\usage{
lrm(formula, data, subset, na.action=na.delete, method="lrm.fit",
    model=FALSE, x=FALSE, y=FALSE, linear.predictors=TRUE, se.fit=FALSE, 
    penalty=0, penalty.matrix, tol=1e-7, 
    strata.penalty=0, var.penalty=c('simple','sandwich'),
    weights, normwt, \dots)

\method{print}{lrm}(x, digits=4, strata.coefs=FALSE,
    coefs=TRUE, latex=FALSE, title='Logistic Regression Model', \dots)
}
\arguments{
\item{formula}{
a formula object. An \code{offset} term can be included. The offset causes
fitting of a model such as \eqn{logit(Y=1) = X\beta + W}, where \eqn{W} is the
offset variable having no estimated coefficient.
The response variable can be any data type; \code{lrm} converts it
in alphabetic or numeric order to an S factor variable and
recodes it 0,1,2,\dots internally. 
}
\item{data}{
data frame to use. Default is the current frame.
}
\item{subset}{
logical expression or vector of subscripts defining a subset of
observations to analyze
}
\item{na.action}{
function to handle \code{NA}s in the data. Default is \code{na.delete}, which
deletes any observation having response or predictor missing, while
preserving the attributes of the predictors and maintaining frequencies
of deletions due to each variable in the model.  
This is usually specified using \code{options(na.action="na.delete")}.
}
\item{method}{
name of fitting function. Only allowable choice at present is \code{lrm.fit}.
}
\item{model}{
causes the model frame to be returned in the fit object
}
\item{x}{
causes the expanded design matrix (with missings excluded)
to be returned under the name \code{x}.  For \code{print}, an object
created by \code{lrm}.
}
\item{y}{
causes the response variable (with missings excluded) to be returned
under the name \code{y}.
}
\item{linear.predictors}{
causes the predicted X beta (with missings excluded) to be returned
under the name \code{linear.predictors}.  When the response variable has
more than two levels, the first intercept is used.
}
\item{se.fit}{
causes the standard errors of the fitted values to be returned under
the name \code{se.fit}.
}
\item{penalty}{
The penalty factor subtracted from the log likelihood is
\eqn{0.5 \beta' P \beta}, where \eqn{\beta} is the vector of regression
coefficients other than intercept(s), and \eqn{P} is
\code{penalty factors * penalty.matrix} and \code{penalty.matrix} is
defined below.  The default is \code{penalty=0} implying that ordinary 
unpenalized maximum likelihood estimation is used.
If \code{penalty} is a scalar, it is assumed to be a penalty factor that
applies 
to all non-intercept parameters in the model.  Alternatively, specify a
list to penalize different types of model terms by differing amounts.
The elements in this list are named \code{simple, nonlinear, interaction} and
\code{nonlinear.interaction}.  If you omit elements on the right of this
series, values are inherited from elements on the left.  Examples:
\code{penalty=list(simple=5, nonlinear=10)} uses a penalty factor of 10
for nonlinear or interaction terms.  
\code{penalty=list(simple=0, nonlinear=2, nonlinear.interaction=4)} does not
penalize linear main effects, uses a penalty factor of 2 for nonlinear or
interaction effects (that are not both), and 4 for nonlinear interaction
effects.
}
\item{penalty.matrix}{
specifies the symmetric penalty matrix for non-intercept terms.
The default matrix for continuous predictors has
the variance of the columns of the design matrix in its diagonal elements
so that the penalty to the log likelhood is unitless.  For main effects
for categorical predictors with \eqn{c} categories, the rows and columns of
the matrix contain a \eqn{c-1 \times c-1} sub-matrix that is used to
compute the 
sum of squares about the mean of the \eqn{c} parameter values (setting the
parameter to zero for the reference cell) as the penalty component
for that predictor.  This makes the penalty independent of the choice of
the reference cell.  If you specify \code{penalty.matrix}, you may set
the rows and columns for certain parameters to zero so as to not
penalize those parameters.
Depending on \code{penalty}, some elements of \code{penalty.matrix} may
be overridden automatically by setting them to zero.
The penalty matrix that is used in the actual fit is 
\eqn{penalty \times diag(pf) \times penalty.matrix \times diag(pf)},
where \eqn{pf} is the vector 
of square roots of penalty factors computed from \code{penalty} by
\code{Penalty.setup} in \code{rmsMisc}.  If you specify \code{penalty.matrix}
you must specify a nonzero value of \code{penalty} or no penalization will be
done.
}
\item{tol}{singularity criterion (see \code{lrm.fit})}
\item{strata.penalty}{scalar penalty factor for the stratification
  factor, for the experimental \code{strat} variable}
\item{var.penalty}{
the type of variance-covariance matrix to be stored in the \code{var}
component of the fit when penalization is used.  The default is the
inverse of the penalized information matrix.  Specify
\code{var.penalty="sandwich"} to use the sandwich estimator (see below
under \code{var}), which limited simulation studies have shown yields
variances estimates that are too low.
}
\item{weights}{
  a vector (same length as \code{y}) of possibly fractional case weights
}
\item{normwt}{
 set to \code{TRUE} to scale \code{weights} so they sum to the length of
 \code{y}; useful for sample surveys as opposed to the default of
 frequency weighting 
 }
\item{\dots}{arguments that are passed to \code{lrm.fit}}
\item{digits}{number of significant digits to use}
\item{strata.coefs}{set to \code{TRUE} to print the (experimental)
  strata coefficients}
\item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
  of model coefficients, standard errors, etc.  Specify \code{coefs=n}
  to print only the first \code{n} regression coefficients in the
  model.}
\item{latex}{a logical value indicating whether information should be
  formatted as plain text or as LaTeX markup}
\item{title}{a character string title to be passed to \code{prModFit}}
}
\value{
The returned fit object of \code{lrm} contains the following components
in addition to the ones mentioned under the optional arguments.

\item{call}{
calling expression
}
\item{freq}{
table of frequencies for \code{Y} in order of increasing \code{Y}
}
\item{stats}{
vector with the following elements: number of observations used in the
fit, maximum absolute value of first
derivative of log likelihood, model likelihood ratio
\eqn{\chi^2}{chi-square}, d.f., 
\eqn{P}-value, \eqn{c} index (area under ROC curve), Somers' \eqn{D_{xy}},
Goodman-Kruskal \eqn{\gamma}{gamma}, Kendall's \eqn{\tau_a}{tau-a} rank
correlations 
between predicted probabilities and observed response, the
Nagelkerke \eqn{R^2} index, the Brier score computed with respect to
\eqn{Y >} its lowest level, the \eqn{g}-index, \eqn{gr} (the
\eqn{g}-index on the odds ratio scale), and \eqn{gp} (the \eqn{g}-index
on the probability scale using the same cutoff used for the Brier
score).  Probabilities are rounded to the nearest 0.002 
in the computations or rank correlation indexes.
In the case of penalized estimation, the \code{"Model L.R."} is computed
without the penalty factor, and \code{"d.f."} is the effective d.f. from
Gray's (1992) Equation 2.9.
The \eqn{P}-value uses this corrected model
L.R. \eqn{\chi^2}{chi-square} and corrected d.f. 
The score chi-square statistic uses first derivatives which contain
penalty components.
}
\item{fail}{
set to \code{TRUE} if convergence failed (and \code{maxiter>1})
}
\item{coefficients}{estimated parameters}
\item{var}{
estimated variance-covariance matrix (inverse of information matrix).
If \code{penalty>0}, \code{var} is either the inverse of the penalized
information matrix (the default, if \code{var.penalty="simple"}) or the
sandwich-type variance - covariance
matrix estimate (Gray Eq. 2.6) if \code{var.penalty="sandwich"}.  For the
latter case the simple information-matrix - based variance
matrix is returned under the name \code{var.from.info.matrix}.
}
\item{effective.df.diagonal}{
is returned if \code{penalty>0}.  It is the vector whose sum is the effective
d.f. of the model (counting intercept terms).
}
\item{u}{
vector of first derivatives of log-likelihood
}
\item{deviance}{
-2 log likelihoods (counting penalty components)
When an offset variable is present, three
deviances are computed: for intercept(s) only, for
intercepts+offset, and for intercepts+offset+predictors.
When there is no offset variable, the vector contains deviances for
the intercept(s)-only model and the model with intercept(s) and predictors.
}
\item{est}{
vector of column numbers of \code{X} fitted (intercepts are not counted)
}
\item{non.slopes}{
number of intercepts in model
}
\item{penalty}{
see above
}
\item{penalty.matrix}{
the penalty matrix actually used in the estimation
}
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Le Cessie S, Van Houwelingen JC: Ridge estimators in logistic regression.
Applied Statistics 41:191--201, 1992.


Verweij PJM, Van Houwelingen JC: Penalized likelihood in Cox regression.
Stat in Med 13:2427--2436, 1994.


Gray RJ: Flexible methods for analyzing survival data using splines,
with applications to breast cancer prognosis.  JASA 87:942--951, 1992.


Shao J: Linear model selection by cross-validation.  JASA 88:486--494, 1993.


Verweij PJM, Van Houwelingen JC: Crossvalidation in survival analysis.
Stat in Med 12:2305--2314, 1993.


Harrell FE: Model uncertainty, penalization, and parsimony.  ISCB
Presentation on UVa Web page, 1998.
}
\seealso{
\code{\link{lrm.fit}}, \code{\link{predict.lrm}},
\code{\link{rms.trans}}, \code{\link{rms}}, \code{\link{glm}},
\code{\link{latex.lrm}}, 
\code{\link{residuals.lrm}}, \code{\link[Hmisc]{na.delete}},
\code{\link[Hmisc]{na.detail.response}},
\code{\link{pentrace}}, \code{\link{rmsMisc}}, \code{\link{vif}},
\code{\link{cr.setup}}, \code{\link{predab.resample}},
\code{\link{validate.lrm}}, \code{\link{calibrate}},
\code{\link{Mean.lrm}}, \code{\link{gIndex}}, \code{\link{prModFit}}
}
\examples{
#Fit a logistic model containing predictors age, blood.pressure, sex
#and cholesterol, with age fitted with a smooth 5-knot restricted cubic 
#spline function and a different shape of the age relationship for males 
#and females.  As an intermediate step, predict mean cholesterol from
#age using a proportional odds ordinal logistic model
#
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

#To use prop. odds model, avoid using a huge number of intercepts by
#grouping cholesterol into 40-tiles
ch <- cut2(cholesterol, g=40, levels.mean=TRUE) # use mean values in intervals
table(ch)
f <- lrm(ch ~ age)
print(f, latex=TRUE, coefs=4)
m <- Mean(f)    # see help file for Mean.lrm
d <- data.frame(age=seq(0,90,by=10))
m(predict(f, d))
# Repeat using ols
f <- ols(cholesterol ~ age)
predict(f, d)

# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
     (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)
cholesterol[1:3] <- NA   # 3 missings, at random

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)
#      x=TRUE, y=TRUE allows use of resid(), which.influence below
#      could define d <- datadist(fit) after lrm(), but data distribution
#      summary would not be stored with fit, so later uses of Predict
#      or summary.rms would require access to the original dataset or
#      d or specifying all variable values to summary, Predict, nomogram
anova(fit)
p <- Predict(fit, age, sex)
plot(p)
plot(Predict(fit, age=20:70, sex="male"))   # need if datadist not used
print(cbind(resid(fit,"dfbetas"), resid(fit,"dffits"))[1:20,])
which.influence(fit, .3)
# latex(fit)                       #print nice statement of fitted model
#
#Repeat this fit using penalized MLE, penalizing complex terms
#(for nonlinear or interaction effects)
#
fitp <- update(fit, penalty=list(simple=0,nonlinear=10), x=TRUE, y=TRUE)
effective.df(fitp)
# or lrm(y ~ \dots, penalty=\dots)


#Get fits for a variety of penalties and assess predictive accuracy 
#in a new data set.  Program efficiently so that complex design 
#matrices are only created once.


set.seed(201)
x1 <- rnorm(500)
x2 <- rnorm(500)
x3 <- sample(0:1,500,rep=TRUE)
L  <- x1+abs(x2)+x3
y  <- ifelse(runif(500)<=plogis(L), 1, 0)
new.data <- data.frame(x1,x2,x3,y)[301:500,]
#
for(penlty in seq(0,.15,by=.005)) {
  if(penlty==0) {
    f <- lrm(y ~ rcs(x1,4)+rcs(x2,6)*x3, subset=1:300, x=TRUE, y=TRUE)
    # True model is linear in x1 and has no interaction
    X <- f$x    # saves time for future runs - don't have to use rcs etc.
    Y <- f$y    # this also deletes rows with NAs (if there were any)
    penalty.matrix <- diag(diag(var(X)))
    Xnew <- predict(f, new.data, type="x", incl.non.slopes=FALSE)  
    # expand design matrix for new data
    Ynew <- new.data$y
  } else f <- lrm.fit(X,Y, penalty.matrix=penlty*penalty.matrix)
#
  cat("\nPenalty :",penlty,"\n")
  pred.logit <- f$coef[1] + (Xnew \%*\% f$coef[-1])
  pred <- plogis(pred.logit)
  C.index <- somers2(pred, Ynew)["C"]
  Brier   <- mean((pred-Ynew)^2)
  Deviance<- -2*sum( Ynew*log(pred) + (1-Ynew)*log(1-pred) )
  cat("ROC area:",format(C.index),"   Brier score:",format(Brier),
      "   -2 Log L:",format(Deviance),"\n")
}
#penalty=0.045 gave lowest -2 Log L, Brier, ROC in test sample for S+
#
#Use bootstrap validation to estimate predictive accuracy of
#logistic models with various penalties
#To see how noisy cross-validation estimates can be, change the
#validate(f, \dots) to validate(f, method="cross", B=10) for example.
#You will see tremendous variation in accuracy with minute changes in
#the penalty.  This comes from the error inherent in using 10-fold
#cross validation but also because we are not fixing the splits.  
#20-fold cross validation was even worse for some
#indexes because of the small test sample size.  Stability would be
#obtained by using the same sample splits for all penalty values 
#(see above), but then we wouldn't be sure that the choice of the 
#best penalty is not specific to how the sample was split.  This
#problem is addressed in the last example.
#
penalties <- seq(0,.7,by=.1)   # really use by=.02
index <- matrix(NA, nrow=length(penalties), ncol=11,
	        dimnames=list(format(penalties),
          c("Dxy","R2","Intercept","Slope","Emax","D","U","Q","B","g","gp")))
i <- 0
for(penlty in penalties)
{
  cat(penlty, "")
  i <- i+1
  if(penlty==0)
    {
    f <- lrm(y ~ rcs(x1,4)+rcs(x2,6)*x3, x=TRUE, y=TRUE)  # fit whole sample
    X <- f$x
    Y <- f$y
    penalty.matrix <- diag(diag(var(X)))   # save time - only do once
    }
  else
   f <- lrm(Y ~ X, penalty=penlty,
            penalty.matrix=penalty.matrix, x=TRUE,y=TRUE)
  val <- validate(f, method="boot", B=20)  # use larger B in practice
  index[i,] <- val[,"index.corrected"]
}
par(mfrow=c(3,3))
for(i in 1:9)
{
  plot(penalties, index[,i], 
       xlab="Penalty", ylab=dimnames(index)[[2]][i])
  lines(lowess(penalties, index[,i]))
}
options(datadist=NULL)

# Example of weighted analysis
x <- 1:5
y <- c(0,1,0,1,0)
reps <- c(1,2,3,2,1)
lrm(y ~ x, weights=reps)
x <- rep(x, reps)
y <- rep(y, reps)
lrm(y ~ x)   # same as above

#
#Study performance of a modified AIC which uses the effective d.f.
#See Verweij and Van Houwelingen (1994) Eq. (6).  Here AIC=chisq-2*df.
#Also try as effective d.f. equation (4) of the previous reference.
#Also study performance of Shao's cross-validation technique (which was
#designed to pick the "right" set of variables, and uses a much smaller
#training sample than most methods).  Compare cross-validated deviance
#vs. penalty to the gold standard accuracy on a 7500 observation dataset.
#Note that if you only want to get AIC or Schwarz Bayesian information
#criterion, all you need is to invoke the pentrace function.
#NOTE: the effective.df( ) function is used in practice
#
\dontrun{
for(seed in c(339,777,22,111,3)){ 
# study performance for several datasets
  set.seed(seed)
  n <- 175; p <- 8
  X <- matrix(rnorm(n*p), ncol=p) # p normal(0,1) predictors
  Coef <- c(-.1,.2,-.3,.4,-.5,.6,-.65,.7)  # true population coefficients
  L <- X \%*\% Coef                 # intercept is zero
  Y <- ifelse(runif(n)<=plogis(L), 1, 0)
  pm <- diag(diag(var(X)))
  #Generate a large validation sample to use as a gold standard
  n.val <- 7500
  X.val <- matrix(rnorm(n.val*p), ncol=p)
  L.val <- X.val \%*\% Coef
  Y.val <- ifelse(runif(n.val)<=plogis(L.val), 1, 0)
  #
  Penalty <- seq(0,30,by=1)
  reps <- length(Penalty)
  effective.df <- effective.df2 <- aic <- aic2 <- deviance.val <- 
    Lpenalty <- single(reps)
  n.t <- round(n^.75)
  ncv <- c(10,20,30,40)     # try various no. of reps in cross-val.
  deviance <- matrix(NA,nrow=reps,ncol=length(ncv))
  #If model were complex, could have started things off by getting X, Y
  #penalty.matrix from an initial lrm fit to save time
  #
  for(i in 1:reps) {
    pen <- Penalty[i]
    cat(format(pen),"")
    f.full <- lrm.fit(X, Y, penalty.matrix=pen*pm)
    Lpenalty[i] <- pen* t(f.full$coef[-1]) \%*\% pm \%*\% f.full$coef[-1]
    f.full.nopenalty <- lrm.fit(X, Y, initial=f.full$coef, maxit=1)
    info.matrix.unpenalized <- solve(f.full.nopenalty$var)
    effective.df[i] <- sum(diag(info.matrix.unpenalized \%*\% f.full$var)) - 1
    lrchisq <- f.full.nopenalty$stats["Model L.R."]
    # lrm does all this penalty adjustment automatically (for var, d.f.,
    # chi-square)
    aic[i] <- lrchisq - 2*effective.df[i]
    #
    pred <- plogis(f.full$linear.predictors)
    score.matrix <- cbind(1,X) * (Y - pred)
    sum.u.uprime <- t(score.matrix) \%*\% score.matrix
    effective.df2[i] <- sum(diag(f.full$var \%*\% sum.u.uprime))
    aic2[i] <- lrchisq - 2*effective.df2[i]
    #
    #Shao suggested averaging 2*n cross-validations, but let's do only 40
    #and stop along the way to see if fewer is OK
    dev <- 0
    for(j in 1:max(ncv)) {
      s    <- sample(1:n, n.t)
      cof  <- lrm.fit(X[s,],Y[s], 
                      penalty.matrix=pen*pm)$coef
      pred <- cof[1] + (X[-s,] \%*\% cof[-1])
      dev <- dev -2*sum(Y[-s]*pred + log(1-plogis(pred)))
      for(k in 1:length(ncv)) if(j==ncv[k]) deviance[i,k] <- dev/j
    }
    #
    pred.val <- f.full$coef[1] + (X.val \%*\% f.full$coef[-1])
    prob.val <- plogis(pred.val)
    deviance.val[i] <- -2*sum(Y.val*pred.val + log(1-prob.val))
  }
  postscript(hor=TRUE)   # along with graphics.off() below, allow plots
  par(mfrow=c(2,4))   # to be printed as they are finished
  plot(Penalty, effective.df, type="l")
  lines(Penalty, effective.df2, lty=2)
  plot(Penalty, Lpenalty, type="l")
  title("Penalty on -2 log L")
  plot(Penalty, aic, type="l")
  lines(Penalty, aic2, lty=2)
  for(k in 1:length(ncv)) {
    plot(Penalty, deviance[,k], ylab="deviance")
    title(paste(ncv[k],"reps"))
    lines(supsmu(Penalty, deviance[,k]))
  }
  plot(Penalty, deviance.val, type="l")
  title("Gold Standard (n=7500)")
  title(sub=format(seed),adj=1,cex=.5)
  graphics.off()
}
}
#The results showed that to obtain a clear picture of the penalty-
#accuracy relationship one needs 30 or 40 reps in the cross-validation.
#For 4 of 5 samples, though, the super smoother was able to detect
#an accurate penalty giving the best (lowest) deviance using 10-fold
#cross-validation.  Cross-validation would have worked better had
#the same splits been used for all penalties.
#The AIC methods worked just as well and are much quicker to compute.
#The first AIC based on the effective d.f. in Gray's Eq. 2.9
#(Verweij and Van Houwelingen (1994) Eq. 5 (note typo)) worked best.
}
\keyword{category}
\keyword{models}
\concept{logistic regression model}
\concept{ordinal logistic model}
\concept{proportional odds model}
\concept{continuation ratio model}
\concept{ordinal response}
