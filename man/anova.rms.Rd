\name{anova.rms}
\alias{anova.rms}
\alias{print.anova.rms}
\alias{plot.anova.rms}
\alias{latex.anova.rms}
\alias{html.anova.rms}
\title{Analysis of Variance (Wald and F Statistics)}
\description{
The \code{anova} function automatically tests most meaningful hypotheses
in a design. For example, suppose that age and cholesterol are
predictors, and that a general interaction is modeled using a restricted
spline surface. \code{anova} prints Wald statistics (\eqn{F} statistics
for an \code{ols} fit) for testing linearity of age, linearity of
cholesterol, age effect (age + age by cholesterol interaction),
cholesterol effect (cholesterol + age by cholesterol interaction),
linearity of the age by cholesterol interaction (i.e., adequacy of the
simple age * cholesterol 1 d.f. product), linearity of the interaction
in age alone, and linearity of the interaction in cholesterol
alone. Joint tests of all interaction terms in the model and all
nonlinear terms in the model are also performed.  For any multiple
d.f. effects for continuous variables that were not modeled through
\code{rcs}, \code{pol}, \code{lsp}, etc., tests of linearity will be
omitted.  This applies to matrix predictors produced by e.g.
\code{poly} or \code{ns}.  \code{print.anova.rms} is the printing
method.  \code{plot.anova.rms} draws dot charts depicting the importance
of variables in the model, as measured by Wald \eqn{\chi^2}{chi-square},
\eqn{\chi^2}{chi-square} minus d.f., AIC, \eqn{P}-values, partial
\eqn{R^2}, \eqn{R^2} for the whole model after deleting the effects in
question, or proportion of overall model \eqn{R^2} that is due to each
predictor.  \code{latex.anova.rms} is the \code{latex} method.  It
substitutes Greek/math symbols in column headings, uses boldface for
\code{TOTAL} lines, and constructs a caption.  Then it passes the result
to \code{latex.default} for conversion to LaTeX.

For Bayesian models such as \code{blrm}, \code{anova} computes relative
explained variation indexes (REV) based on approximate Wald statistics.
This uses the variance-covariance matrix of all of the posterior draws,
and the individual draws of betas, plus an overall summary from the
posterior mode/mean/median beta.  Wald chi-squares assuming multivariate
normality of betas are computed just as with frequentist models, and for
each draw (or for the summary) the ratio of the partial Wald chi-square
to the total Wald statistic for the model is computed as REV.

The \code{print} method calls \code{latex} or \code{html} methods
depending on \code{options(prType=)}.  For
\code{latex} a \code{table} environment is not used and an ordinary
\code{tabular} is produced.  When using html with Quarto or RMarkdown,
  \code{results='asis'} need not be written in the chunk header.

\code{html.anova.rms} just calls \code{latex.anova.rms}.
}
\usage{
\method{anova}{rms}(object, \ldots, main.effect=FALSE, tol=1e-9,
      test=c('F','Chisq'), india=TRUE, indnl=TRUE, ss=TRUE,
      vnames=c('names','labels'),
      posterior.summary=c('mean', 'median', 'mode'), ns=500, cint=0.95)

\method{print}{anova.rms}(x,
      which=c('none','subscripts','names','dots'),
      table.env=FALSE, \dots)

\method{plot}{anova.rms}(x,
     what=c("chisqminusdf","chisq","aic","P","partial R2","remaining R2",
            "proportion R2", "proportion chisq"),
     xlab=NULL, pch=16,
     rm.totals=TRUE, rm.ia=FALSE, rm.other=NULL, newnames,
     sort=c("descending","ascending","none"), margin=c('chisq','P'),
     pl=TRUE, trans=NULL, ntrans=40, height=NULL, width=NULL, \dots)

\method{latex}{anova.rms}(object, title, dec.chisq=2,
      dec.F=2, dec.ss=NA, dec.ms=NA, dec.P=4, dec.REV=3, table.env=TRUE,
      caption=NULL, fontsize=1, params, \dots)

\method{html}{anova.rms}(object, \dots)
}
\arguments{
\item{object}{
a \code{rms} fit object.  \code{object} must
allow \code{vcov} to return the variance-covariance matrix.  For
\code{latex} is the result of \code{anova}.
}
\item{\dots}{
If omitted, all variables are tested, yielding tests for individual factors
and for pooled effects. Specify a subset of the variables to obtain tests
for only those factors, with a pooled Wald tests for the combined effects
of all factors listed. Names may be abbreviated.  For example, specify
\code{anova(fit,age,cholesterol)} to get a Wald statistic for testing the joint
importance of age, cholesterol, and any factor interacting with them.

Can be optional graphical parameters to send to
\code{dotchart2}, or other parameters to send to \code{latex.default}.
Ignored for \code{print}.

For \code{html.anova.rms} the arguments are passed to \code{latex.anova.rms}.
}
\item{main.effect}{
Set to \code{TRUE} to print the (usually meaningless) main effect tests even when
the factor is involved in an interaction. The default is \code{FALSE}, to print only
the effect of the main effect combined with all interactions involving that
factor.
}
\item{tol}{
singularity criterion for use in matrix inversion
}
\item{test}{
For an \code{ols} fit, set \code{test="Chisq"} to use Wald \eqn{\chi^2} tests rather than F-tests.
}
\item{india}{set to \code{FALSE} to exclude individual tests of
  interaction from the table}
\item{indnl}{set to \code{FALSE} to exclude individual tests of
  nonlinearity from the table}
\item{ss}{
For an \code{ols} fit, set \code{ss=FALSE} to suppress printing partial
sums of squares, mean squares, and the Error SS and MS.
}
\item{vnames}{set to \code{'labels'} to use variable labels rather than
  variable names in the output}
\item{posterior.summary}{specifies whether the posterior mode/mean/median
	beta are to be used as a measure of central tendence of the posterior
	distribution, for use in relative explained variation from Bayesian
	models}
\item{ns}{number of random samples from the posterior draws to use for
	REV highest posterior density intervals}
\item{cint}{HPD interval probability}
\item{x}{for \code{print,plot,text} is the result of \code{anova}.
}
\item{which}{
If \code{which} is not \code{"none"} (the default), \code{print.anova.rms} will
add to the rightmost column of the output the list of parameters being
tested by the hypothesis being tested in the current row.  Specifying
\code{which="subscripts"} causes the subscripts of the regression
coefficients being tested to be printed (with a subscript of one for
the first non-intercept term).  \code{which="names"} prints the names of
the terms being tested, and \code{which="dots"} prints dots for terms being
tested and blanks for those just being adjusted for.
}
\item{what}{
what type of statistic to plot.  The default is the Wald
\eqn{\chi^2}{chi-square}
statistic for each factor (adding in the effect of higher-ordered
factors containing that factor) minus its degrees of freedom.  The
R2 choices for \code{what} only apply to \code{ols} models.
}
\item{xlab}{
x-axis label, default is constructed according to \code{what}.
\code{plotmath} symbols are used for \R, by default.
}
\item{pch}{
character for plotting dots in dot charts.  Default is 16 (solid dot).
}
\item{rm.totals}{
set to \code{FALSE} to keep total \eqn{\chi^2}{chi-square}s (overall, nonlinear, interaction totals)
in the chart.
}
\item{rm.ia}{
set to \code{TRUE} to omit any effect that has \code{"*"} in its name
}
\item{rm.other}{
a list of other predictor names to omit from the chart
}
\item{newnames}{
a list of substitute predictor names to use, after omitting any.
}
\item{sort}{default is to sort bars in descending order of the summary statistic. 
Available options: 'ascending', 'descending', 'none'. 
}
\item{margin}{set to a vector of character strings to write text for
	selected statistics in the right margin of the dot chart.  The
	character strings can be any combination of \code{"chisq"},
	\code{"d.f."}, \code{"P"}, \code{"partial R2"},
	\code{"proportion R2"}, and \code{"proportion chisq"}.
	Default is to not draw any statistics in the margin.  When
	\code{plotly} is in effect, margin values are instead displayed as
	hover text.}
\item{pl}{
set to \code{FALSE} to suppress plotting.  This is useful when you only wish to
analyze the vector of statistics returned.
}
\item{trans}{
  set to a function to apply that transformation to the statistics
  being plotted, and to truncate negative values at zero.  A good choice
  is \code{trans=sqrt}.
}
\item{ntrans}{\code{n} argument to \code{\link{pretty}}, specifying the
	number of values for which to place tick marks.  This should be larger
	than usual because of nonlinear scaling, to provide a sufficient
	number of tick marks on the left (stretched) part of the chi-square
	scale.
}
\item{height,width}{height and width of \code{plotly} plots drawn using
	\code{dotchartp}, in pixels.  Ignored for ordinary plots.  Defaults to
minimum of 400 and 100 + 25 times the number of test statistics displayed.}
\item{title}{
title to pass to \code{latex}, default is name of fit object passed to
\code{anova} prefixed with \code{"anova."}.  For Windows, the default is
\code{"ano"} followed by the first 5 letters of the name of the fit
object.
}
\item{dec.chisq}{
number of places to the right of the decimal place for typesetting
\eqn{\chi^2}{chi-square} values (default is \code{2}).  Use zero for integer, \code{NA} for
floating point.
}
\item{dec.F}{
digits to the right for \eqn{F} statistics (default is \code{2})
}
\item{dec.ss}{
digits to the right for sums of squares (default is \code{NA}, indicating
floating point)
}
\item{dec.ms}{
digits to the right for mean squares (default is \code{NA})
}
\item{dec.P}{digits to the right for \eqn{P}-values}
\item{dec.REV}{digits to the right for REV}
\item{table.env}{see \code{\link[Hmisc]{latex}}}
\item{caption}{caption for table if \code{table.env} is \code{TRUE}.
  Default is constructed from the response variable.}
\item{fontsize}{font size for html output; default is 1 for \code{1em}}
\item{params}{used internally when called through print.}
}
\value{
\code{anova.rms} returns a matrix of class \code{anova.rms} containing factors
as rows and \eqn{\chi^2}{chi-square}, d.f., and \eqn{P}-values as
columns (or d.f., partial \eqn{SS, MS, F, P}).  An attribute
\code{vinfo} provides list of variables involved in each row and the
type of test done.
\code{plot.anova.rms} invisibly returns the vector of quantities
plotted.  This vector has a names attribute describing the terms for
which the statistics in the vector are calculated.
}
\details{
If the statistics being plotted with \code{plot.anova.rms} are few in
number and one of them is negative or zero, \code{plot.anova.rms}
will quit because of an error in \code{dotchart2}.

The \code{latex} method requires LaTeX packages \code{relsize} and
\code{needspace}.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\section{Side Effects}{
\code{print} prints, \code{latex} creates a
file with a name of the form \code{"title.tex"} (see the \code{title} argument above).
}
\seealso{
\code{\link{rms}}, \code{\link{rmsMisc}}, \code{\link{lrtest}},
\code{\link{rms.trans}}, \code{\link{summary.rms}}, \code{\link{plot.Predict}},
\code{\link{ggplot.Predict}}, \code{\link[Hmisc]{solvet}},
\code{\link{locator}},
\code{\link[Hmisc]{dotchart2}}, \code{\link[Hmisc]{latex}},
\code{\link[Hmisc]{xYplot}}, \code{\link{anova.lm}},
\code{\link{contrast.rms}}, \code{\link{pantext}}
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
treat <- factor(sample(c('a','b','c'), n,TRUE))
num.diseases <- sample(0:4, n,TRUE)
age <- rnorm(n, 50, 10)
cholesterol <- rnorm(n, 200, 25)
weight <- rnorm(n, 150, 20)
sex <- factor(sample(c('female','male'), n,TRUE))
label(age) <- 'Age'      # label is in Hmisc
label(num.diseases) <- 'Number of Comorbid Diseases'
label(cholesterol) <- 'Total Cholesterol'
label(weight) <- 'Weight, lbs.'
label(sex) <- 'Sex'
units(cholesterol) <- 'mg/dl'   # uses units.default in Hmisc


# Specify population model for log odds that Y=1
L <- .1*(num.diseases-2) + .045*(age-50) +
     (log(cholesterol - 10)-5.2)*(-2*(treat=='a') +
     3.5*(treat=='b')+2*(treat=='c'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)


fit <- lrm(y ~ treat + scored(num.diseases) + rcs(age) +
               log(cholesterol+10) + treat:log(cholesterol+10))
a <- anova(fit)                       # Test all factors
b <- anova(fit, treat, cholesterol)   # Test these 2 by themselves
                                      # to get their pooled effects
a
b
# Add a new line to the plot with combined effects
s <- rbind(a, 'treat+cholesterol'=b['TOTAL',])
class(s) <- 'anova.rms'
plot(s, margin=c('chisq', 'proportion chisq'))

g <- lrm(y ~ treat*rcs(age))
dd <- datadist(treat, num.diseases, age, cholesterol)
options(datadist='dd')
p <- Predict(g, age, treat="b")
s <- anova(g)
tx <- paste(capture.output(s), collapse='\n')
ggplot(p) + annotate('text', x=27, y=3.2, family='mono', label=tx,
                      hjust=0, vjust=1, size=1.5)

plot(s, margin=c('chisq', 'proportion chisq'))
# new plot - dot chart of chisq-d.f. with 2 other stats in right margin
# latex(s)                       # nice printout - creates anova.g.tex
options(datadist=NULL)


# Simulate data with from a given model, and display exactly which
# hypotheses are being tested


set.seed(123)
age <- rnorm(500, 50, 15)
treat <- factor(sample(c('a','b','c'), 500, TRUE))
bp  <- rnorm(500, 120, 10)
y   <- ifelse(treat=='a', (age-50)*.05, abs(age-50)*.08) + 3*(treat=='c') +
       pmax(bp, 100)*.09 + rnorm(500)
f   <- ols(y ~ treat*lsp(age,50) + rcs(bp,4))
print(names(coef(f)), quote=FALSE)
specs(f)
anova(f)
an <- anova(f)
options(digits=3)
print(an, 'subscripts')
print(an, 'dots')


an <- anova(f, test='Chisq', ss=FALSE)
# plot(0:1)                        # make some plot
# tab <- pantext(an, 1.2, .6, lattice=FALSE, fontfamily='Helvetica')
# create function to write table; usually omit fontfamily
# tab()                            # execute it; could do tab(cex=.65)
plot(an)                         # new plot - dot chart of chisq-d.f.
# Specify plot(an, trans=sqrt) to use a square root scale for this plot
# latex(an)                      # nice printout - creates anova.f.tex


## Example to save partial R^2 for all predictors, along with overall
## R^2, from two separate fits, and to combine them with ggplot2

require(ggplot2)
set.seed(1)
n <- 100
x1 <- runif(n)
x2 <- runif(n)
y  <- (x1-.5)^2 + x2 + runif(n)
group <- c(rep('a', n/2), rep('b', n/2))
A <- NULL
for(g in c('a','b')) {
    f <- ols(y ~ pol(x1,2) + pol(x2,2) + pol(x1,2) \%ia\% pol(x2,2),
             subset=group==g)
    a <- plot(anova(f),
              what='partial R2', pl=FALSE, rm.totals=FALSE, sort='none')
    a <- a[-grep('NONLINEAR', names(a))]
    d <- data.frame(group=g, Variable=factor(names(a), names(a)),
                    partialR2=unname(a))
    A <- rbind(A, d)
  }
ggplot(A, aes(x=partialR2, y=Variable)) + geom_point() +
       facet_wrap(~ group) + xlab(ex <- expression(partial~R^2)) +
       scale_y_discrete(limits=rev)
ggplot(A, aes(x=partialR2, y=Variable, color=group)) + geom_point() +
       xlab(ex <- expression(partial~R^2)) +
       scale_y_discrete(limits=rev)

# Suppose that a researcher wants to make a big deal about a variable
# because it has the highest adjusted chi-square.  We use the
# bootstrap to derive 0.95 confidence intervals for the ranks of all
# the effects in the model.  We use the plot method for anova, with
# pl=FALSE to suppress actual plotting of chi-square - d.f. for each
# bootstrap repetition.
# It is important to tell plot.anova.rms not to sort the results, or
# every bootstrap replication would have ranks of 1,2,3,... for the stats.

n <- 300
set.seed(1)
d <- data.frame(x1=runif(n), x2=runif(n),  x3=runif(n),
   x4=runif(n), x5=runif(n), x6=runif(n),  x7=runif(n),
   x8=runif(n), x9=runif(n), x10=runif(n), x11=runif(n),
   x12=runif(n))
d$y <- with(d, 1*x1 + 2*x2 + 3*x3 +  4*x4  + 5*x5 + 6*x6 +
               7*x7 + 8*x8 + 9*x9 + 10*x10 + 11*x11 +
              12*x12 + 9*rnorm(n))

f <- ols(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12, data=d)
B <- 20   # actually use B=1000
ranks <- matrix(NA, nrow=B, ncol=12)
rankvars <- function(fit)
  rank(plot(anova(fit), sort='none', pl=FALSE))
Rank <- rankvars(f)
for(i in 1:B) {
  j <- sample(1:n, n, TRUE)
  bootfit <- update(f, data=d, subset=j)
  ranks[i,] <- rankvars(bootfit)
  }
lim <- t(apply(ranks, 2, quantile, probs=c(.025,.975)))
predictor <- factor(names(Rank), names(Rank))
w <- data.frame(predictor, Rank, lower=lim[,1], upper=lim[,2])
ggplot(w, aes(x=predictor, y=Rank)) + geom_point() + coord_flip() +
  scale_y_continuous(breaks=1:12) +
  geom_errorbar(aes(ymin=lim[,1], ymax=lim[,2]), width=0)
}
\keyword{models}
\keyword{regression}
\keyword{htest}
\keyword{aplot}
\concept{bootstrap}
