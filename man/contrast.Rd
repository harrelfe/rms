\name{contrast.rms}
\alias{contrast}
\alias{contrast.rms}
\alias{print.contrast.rms}
\title{General Contrasts of Regression Coefficients}

\description{
This function computes one or more contrasts of the estimated
regression coefficients in a fit from one of the functions in rms,
along with standard errors, confidence limits, t or Z statistics, P-values.
General contrasts are handled by obtaining the design matrix for two
sets of predictor settings (\code{a}, \code{b}) and subtracting the
corresponding rows of the two design matrics to obtain a new contrast
design matrix for testing the \code{a} - \code{b} differences.  This allows for
quite general contrasts (e.g., estimated differences in means between
a 30 year old female and a 40 year old male).
This can also be used
to obtain a series of contrasts in the presence of interactions (e.g.,
female:male log odds ratios for several ages when the model contains
age by sex interaction).  Another use of \code{contrast} is to obtain
center-weighted (Type III test) and subject-weighted (Type II test)
estimates in a model containing treatment by center interactions.  For
the latter case, you can specify \code{type="average"} and an optional
\code{weights} vector to average the within-center treatment contrasts.
The design contrast matrix computed by \code{contrast.rms} can be used
by other functions.

If \code{usebootcoef=TRUE}, the fit was run through \code{bootcov}, and
\code{conf.type="individual"}, the confidence intervals are bootstrap
nonparametric percentile confidence intervals, basic bootstrap, or BCa
intervals, obtained on contrasts evaluated on all bootstrap samples.

By omitting the \code{b} argument, \code{contrast} can be used to obtain
an average or weighted average of a series of predicted values, along
with a confidence interval for this average.  This can be useful for
"unconditioning" on one of the predictors (see the next to last
example).

Specifying \code{type="joint"}, and specifying at least as many contrasts
as needed to span the space of a complex test, one can make
multiple degree of freedom tests flexibly and simply.  Redundant
contrasts will be ignored in the joint test.  See the examples below.
These include an example of an "incomplete interaction test" involving
only two of three levels of a categorical variable (the test also tests
the main effect).

When more than one contrast is computed, the list created by
\code{contrast.rms} is suitable for plotting (with error bars or bands)
with \code{xYplot} or \code{Dotplot} (see the last example before the
\code{type="joint"} examples).
}
\usage{
contrast(fit, \dots)
\method{contrast}{rms}(fit, a, b, cnames=NULL, 
         type=c("individual", "average", "joint"),
         conf.type=c("individual","simultaneous"), usebootcoef=TRUE,
         boot.type=c("percentile","bca","basic"),
         weights="equal", conf.int=0.95, tol=1e-7, expand=TRUE, \dots)

\method{print}{contrast.rms}(x, X=FALSE, fun=function(u)u, jointonly=FALSE, ...)
}
\arguments{
\item{fit}{
a fit of class \code{"rms"}
}
\item{a}{
a list containing settings for all predictors that you do not wish to
set to default (adjust-to) values.  Usually you will specify two
variables in this list, one set to a constant and one to a sequence of
values, to obtain contrasts for the sequence of values of an
interacting factor.  The \code{gendata} function will generate the
necessary combinations and default values for unspecified predictors,
depending on the \code{expand} argument.
}
\item{b}{
another list that generates the same number of observations as \code{a},
unless one of the two lists generates only one observation.  In that
case, the design matrix generated from the shorter list will have its
rows replicated so that the contrasts assess several differences
against the one set of predictor values.  This is useful for comparing
multiple treatments with control, for example.  If \code{b} is missing, the
design matrix generated from \code{a} is analyzed alone.
}
\item{cnames}{
vector of character strings naming the contrasts when
\code{type!="average"}.  Usually \code{cnames} is not necessary as
\code{contrast.rms} tries to name the contrasts by examining which
predictors are varying consistently in the two lists.  \code{cnames} will
be needed when you contrast "non-comparable" settings, e.g., you
compare \code{list(treat="drug", age=c(20,30))} with
\code{list(treat="placebo"), age=c(40,50))}
}
\item{type}{
set \code{type="average"} to average the individual contrasts (e.g., to
obtain a Type II or III contrast).  Set \code{type="joint"} to jointly
test all non-redundant contrasts with a multiple degree of freedom test
and no averaging.
}
\item{conf.type}{
  The default type of confidence interval computed for a given
  individual (1 d.f.) contrast is a pointwise confidence interval.  Set
  \code{conf.type="simultaneous"} to use the \code{multcomp} package's
  \code{glht} and \code{confint} functions to compute confidence
  intervals with simultaneous (family-wise) coverage, thus adjusting for
  multiple comparisons.  Note that individual P-values are not adjusted
  for multiplicity.
}
\item{usebootcoef}{
  If \code{fit} was the result of \code{bootcov} but you want to use the
  bootstrap covariance matrix instead of the nonparametric percentile,
  basic, or BCa  method for confidence intervals (which uses all the bootstrap
  coefficients), specify \code{usebootcoef=FALSE}.}
\item{boot.type}{set to \code{'bca'} to compute BCa confidence
	limits or \code{'basic'} to use the basic bootstrap.  The default is
	to compute percentile intervals} 
\item{weights}{
a numeric vector, used when \code{type="average"}, to obtain weighted contrasts
}
\item{conf.int}{
confidence level for confidence intervals for the contrasts
}
\item{tol}{tolerance for \code{qr} function for determining which
  contrasts are redundant, and for inverting the covariance matrix
  involved in a joint test}
\item{expand}{set to \code{FALSE} to have \code{gendata} not generate
  all possible combinations of predictor settings.  This is useful when
  getting contrasts over irregular predictor settings.}
\item{\dots}{unused}
\item{x}{result of \code{contrast}}
\item{X}{
set \code{X=TRUE} to  print design matrix used in computing the contrasts (or
the average contrast)
}
\item{fun}{
a function to transform the contrast, SE, and lower and upper
confidence limits before printing.  For example, specify \code{fun=exp} to
anti-log them for logistic models.
}
\item{jointonly}{set to \code{FALSE} to omit printing of individual
  contrasts}
}
\value{
a list of class \code{"contrast.rms"} containing the elements
\code{Contrast}, \code{SE}, \code{Z}, \code{var}, \code{df.residual}
\code{Lower}, \code{Upper}, \code{Pvalue}, \code{X}, \code{cnames}, \code{redundant}, which denote the contrast
estimates, standard errors, Z or t-statistics, variance matrix,
residual degrees of freedom (this is \code{NULL} if the model was not
\code{ols}), lower and upper confidence limits, 2-sided P-value, design
matrix, contrast names (or \code{NULL}), and a logical vector denoting
which contrasts are redundant with the other contrasts.  If there are
any redundant contrasts, when the results of \code{contrast} are
printed, and asterisk is printed at the start of the corresponding
lines.  The object also contains \code{ctype} indicating what method was
used for compute confidence intervals.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University School of Medicine\cr
f.harrell@vanderbilt.edu
}
\seealso{
  \code{\link{Predict}}, \code{\link{gendata}}, \code{\link{bootcov}},
  \code{\link{summary.rms}}, \code{\link{anova.rms}},
}
\examples{
set.seed(1)
age <- rnorm(200,40,12)
sex <- factor(sample(c('female','male'),200,TRUE))
logit <- (sex=='male') + (age-40)/5
y <- ifelse(runif(200) <= plogis(logit), 1, 0)
f <- lrm(y ~ pol(age,2)*sex)
# Compare a 30 year old female to a 40 year old male
# (with or without age x sex interaction in the model)
contrast(f, list(sex='female', age=30), list(sex='male', age=40))


# For a model containing two treatments, centers, and treatment
# x center interaction, get 0.95 confidence intervals separately
# by cente
center <- factor(sample(letters[1:8],500,TRUE))
treat  <- factor(sample(c('a','b'),  500,TRUE))
y      <- 8*(treat=='b') + rnorm(500,100,20)
f <- ols(y ~ treat*center)


lc <- levels(center)
contrast(f, list(treat='b', center=lc),
            list(treat='a', center=lc))


# Get 'Type III' contrast: average b - a treatment effect over
# centers, weighting centers equally (which is almost always
# an unreasonable thing to do)
contrast(f, list(treat='b', center=lc),
            list(treat='a', center=lc),
         type='average')


# Get 'Type II' contrast, weighting centers by the number of
# subjects per center.  Print the design contrast matrix used.
k <- contrast(f, list(treat='b', center=lc),
                 list(treat='a', center=lc),
              type='average', weights=table(center))
print(k, X=TRUE)
# Note: If other variables had interacted with either treat 
# or center, we may want to list settings for these variables
# inside the list()'s, so as to not use default settings


# For a 4-treatment study, get all comparisons with treatment 'a'
treat  <- factor(sample(c('a','b','c','d'),  500,TRUE))
y      <- 8*(treat=='b') + rnorm(500,100,20)
dd     <- datadist(treat,center); options(datadist='dd')
f <- ols(y ~ treat*center)
lt <- levels(treat)
contrast(f, list(treat=lt[-1]),
            list(treat=lt[ 1]),
         cnames=paste(lt[-1],lt[1],sep=':'), conf.int=1-.05/3)


# Compare each treatment with average of all others
for(i in 1:length(lt)) {
  cat('Comparing with',lt[i],'\n\n')
  print(contrast(f, list(treat=lt[-i]),
                    list(treat=lt[ i]), type='average'))
}
options(datadist=NULL)

# Six ways to get the same thing, for a variable that
# appears linearly in a model and does not interact with
# any other variables.  We estimate the change in y per
# unit change in a predictor x1.  Methods 4, 5 also
# provide confidence limits.  Method 6 computes nonparametric
# bootstrap confidence limits.  Methods 2-6 can work
# for models that are nonlinear or non-additive in x1.
# For that case more care is needed in choice of settings
# for x1 and the variables that interact with x1.


\dontrun{
coef(fit)['x1']                            # method 1
diff(predict(fit, gendata(x1=c(0,1))))     # method 2
g <- Function(fit)                         # method 3
g(x1=1) - g(x1=0)
summary(fit, x1=c(0,1))                    # method 4
k <- contrast(fit, list(x1=1), list(x1=0)) # method 5
print(k, X=TRUE)
fit <- update(fit, x=TRUE, y=TRUE)         # method 6
b <- bootcov(fit, B=500)
contrast(fit, list(x1=1), list(x1=0))


# In a model containing age, race, and sex,
# compute an estimate of the mean response for a
# 50 year old male, averaged over the races using
# observed frequencies for the races as weights


f <- ols(y ~ age + race + sex)
contrast(f, list(age=50, sex='male', race=levels(race)),
         type='average', weights=table(race))
}


# Plot the treatment effect (drug - placebo) as a function of age
# and sex in a model in which age nonlinearly interacts with treatment
# for females only

set.seed(1)
n <- 800
treat <- factor(sample(c('drug','placebo'), n,TRUE))
sex   <- factor(sample(c('female','male'),  n,TRUE))
age   <- rnorm(n, 50, 10)
y     <- .05*age + (sex=='female')*(treat=='drug')*.05*abs(age-50) + rnorm(n)
f     <- ols(y ~ rcs(age,4)*treat*sex)
d     <- datadist(age, treat, sex); options(datadist='d')

# show separate estimates by treatment and sex

plot(Predict(f, age, treat, sex='female'))
plot(Predict(f, age, treat, sex='male'))
ages  <- seq(35,65,by=5); sexes <- c('female','male')
w     <- contrast(f, list(treat='drug',    age=ages, sex=sexes),
                     list(treat='placebo', age=ages, sex=sexes))
# add conf.type="simultaneous" to adjust for having done 14 contrasts
xYplot(Cbind(Contrast, Lower, Upper) ~ age | sex, data=w,
       ylab='Drug - Placebo')
xYplot(Cbind(Contrast, Lower, Upper) ~ age, groups=sex, data=w,
       ylab='Drug - Placebo', method='alt bars')
options(datadist=NULL)


# Examples of type='joint' contrast tests

set.seed(1)
x1 <- rnorm(100)
x2 <- factor(sample(c('a','b','c'), 100, TRUE))
dd <- datadist(x1, x2); options(datadist='dd')
y  <- x1 + (x2=='b') + rnorm(100)

# First replicate a test statistic from anova()

f <- ols(y ~ x2)
anova(f)
contrast(f, list(x2=c('b','c')), list(x2='a'), type='joint')

# Repeat with a redundancy; compare a vs b, a vs c, b vs c

contrast(f, list(x2=c('a','a','b')), list(x2=c('b','c','c')), type='joint')

# Get a test of association of a continuous predictor with y
# First assume linearity, then cubic

f <- lrm(y>0 ~ x1 + x2)
anova(f)
contrast(f, list(x1=1), list(x1=0), type='joint')  # a minimum set of contrasts
xs <- seq(-2, 2, length=20)
contrast(f, list(x1=0), list(x1=xs), type='joint')

# All contrasts were redundant except for the first, because of
# linearity assumption

f <- lrm(y>0 ~ pol(x1,3) + x2)
anova(f)
contrast(f, list(x1=0), list(x1=xs), type='joint')
print(contrast(f, list(x1=0), list(x1=xs), type='joint'), jointonly=TRUE)

# All contrasts were redundant except for the first 3, because of
# cubic regression assumption

# Now do something that is difficult to do without cryptic contrast
# matrix operations: Allow each of the three x2 groups to have a different
# shape for the x1 effect where x1 is quadratic.  Test whether there is
# a difference in mean levels of y for x2='b' vs. 'c' or whether
# the shape or slope of x1 is different between x2='b' and x2='c' regardless
# of how they differ when x2='a'.  In other words, test whether the mean
# response differs between group b and c at any value of x1.
# This is a 3 d.f. test (intercept, linear, quadratic effects) and is
# a better approach than subsetting the data to remove x2='a' then
# fitting a simpler model, as it uses a better estimate of sigma from
# all the data.

f <- ols(y ~ pol(x1,2) * x2)
anova(f)
contrast(f, list(x1=xs, x2='b'),
            list(x1=xs, x2='c'), type='joint')

# Note: If using a spline fit, there should be at least one value of
# x1 between any two knots and beyond the outer knots.
options(datadist=NULL)
}
\keyword{htest}
\keyword{models}
\keyword{regression}
