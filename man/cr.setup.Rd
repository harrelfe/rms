\name{cr.setup}
\alias{cr.setup}
\title{Continuation Ratio Ordinal Logistic Setup}
\description{
  Creates several new variables which help set up a dataset with an
  ordinal response variable \eqn{y} for use in fitting a forward continuation
  ratio (CR) model.  The CR model can be fitted with binary logistic
  regression if each input observation is replicated the proper
  number of times according to the \eqn{y} value, a new binary \eqn{y}
  is computed that has at most one \eqn{y=1} per subject,
  and if a \code{cohort} variable is used to define the current
  qualifying condition for a cohort of subjects, e.g., \eqn{y\geq 2}.
  \code{cr.setup} creates the needed auxilliary variables. See
  \code{predab.resample} and \code{validate.lrm} for information about
  validating CR models (e.g., using the bootstrap to sample with
  replacement from the original subjects instead of the records used in
  the fit, validating the model separately for user-specified values of
  \code{cohort}).
}
\usage{
cr.setup(y)
}
\arguments{
\item{y}{
a character, numeric, \code{category}, or \code{factor} vector containing values of
the response variable.  For \code{category} or \code{factor} variables, the
\code{levels} of the variable are assumed to be listed in an ordinal way.
}}
\value{
a list with components \code{y, cohort, subs, reps}.  \code{y} is a new binary
variable that is to be used in the binary logistic fit.  \code{cohort} is 
a \code{factor} vector specifying which cohort condition currently applies.
\code{subs} is a vector of subscripts that can be used to replicate other
variables the same way \code{y} was replicated.  \code{reps} specifies how many
times each original observation was replicated.  \code{y, cohort, subs} are
all the same length and are longer than the original \code{y} vector.
\code{reps} is the same length as the original \code{y} vector.
The \code{subs} vector is suitable for passing to \code{validate.lrm} or \code{calibrate},
which pass this vector under the name \code{cluster} on to \code{predab.resample} so that bootstrapping can be
done by sampling with replacement from the original subjects rather than
from the individual records created by \code{cr.setup}.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
fh@fharrell.com
}
\references{
Berridge DM, Whitehead J: Analysis of failure time data with ordinal
categories of response.  Stat in Med 10:1703--1710, 1991.
}
\seealso{
\code{\link{lrm}}, \code{\link{glm}}, \code{\link{predab.resample}}
}
\examples{
y <- c(NA, 10, 21, 32, 32)
cr.setup(y)


set.seed(171)
y <- sample(0:2, 100, rep=TRUE)
sex <- sample(c("f","m"),100,rep=TRUE)
sex <- factor(sex)
table(sex, y)
options(digits=5)
tapply(y==0, sex, mean)
tapply(y==1, sex, mean)
tapply(y==2, sex, mean)
cohort <- y>=1
tapply(y[cohort]==1, sex[cohort], mean)

u <- cr.setup(y)
Y <- u$y
cohort <- u$cohort
sex <- sex[u$subs]

lrm(Y ~ cohort + sex)
 
f <- lrm(Y ~ cohort*sex)   # saturated model - has to fit all data cells
f

#Prob(y=0|female):
# plogis(-.50078)
#Prob(y=0|male):
# plogis(-.50078+.11301)
#Prob(y=1|y>=1, female):
plogis(-.50078+.31845)
#Prob(y=1|y>=1, male):
plogis(-.50078+.31845+.11301-.07379)

combinations <- expand.grid(cohort=levels(cohort), sex=levels(sex))
combinations
p <- predict(f, combinations, type="fitted")
p
p0 <- p[c(1,3)]
p1 <- p[c(2,4)]
p1.unconditional <- (1 - p0) *p1
p1.unconditional
p2.unconditional <- 1 - p0 - p1.unconditional
p2.unconditional


\dontrun{
dd <- datadist(inputdata)   # do this on non-replicated data
options(datadist='dd')
pain.severity <- inputdata$pain.severity
u <- cr.setup(pain.severity)
# inputdata frame has age, sex with pain.severity
attach(inputdata[u$subs,])  # replicate age, sex
# If age, sex already available, could do age <- age[u$subs] etc., or
# age <- rep(age, u$reps), etc.
y      <- u$y
cohort <- u$cohort
dd     <- datadist(dd, cohort)       # add to dd
f <- lrm(y ~ cohort + age*sex)       # ordinary cont. ratio model
g <- lrm(y ~ cohort*sex + age, x=TRUE,y=TRUE) # allow unequal slopes for
                                     # sex across cutoffs
cal <- calibrate(g, cluster=u$subs, subset=cohort=='all')  
# subs makes bootstrap sample the correct units, subset causes
# Predicted Prob(pain.severity=0) to be checked for calibration
}
}
\keyword{category}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
\concept{continuation ratio model}
\concept{ordinal logistic model}
\concept{ordinal response}
