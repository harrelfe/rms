######################
# Detailed Example 1 #
######################
# May want to first invoke the Hmisc store function
# so that new variables will go into a temporary directory
set.seed(17)  # So can repeat random number sequence
n <- 500

sex    <- factor(sample(c('female','male'), n, rep=TRUE))
age    <- rnorm(n, 50, 10)
sys.bp <- rnorm(n, 120, 7)

# Use two population models, one with a systolic
# blood pressure effect and one without

L    <- ifelse(sex=='female', .1*(pmin(age,50)-50), .005*(age-50)^2)
L.bp <- L + .4*(pmax(sys.bp,120)-120)

dz    <- ifelse(runif(n) <= plogis(L),    1, 0)
dz.bp <- ifelse(runif(n) <= plogis(L.bp), 1, 0)

# Use summary.formula in the Hmisc package to summarize the
# data one predictor at a time

s <- summary(dz.bp ~ age + sex + sys.bp) 
options(digits=3)
print(s)
plot(s)

plsmo(age, dz, group=sex, fun=qlogis, ylim=c(-3,3))
plsmo(age, L,  group=sex, method='raw', add=TRUE, prefix='True', trim=0)
title('Lowess-smoothed Estimates with True Regression Functions')

dd <- datadist(age, sex, sys.bp)
options(datadist='dd')
# can also do: dd <- datadist(dd, newvar)

f <- lrm(dz ~ rcs(age,5)*sex, x=TRUE, y=TRUE)
f
# x=TRUE, y=TRUE for pentrace

fpred <- Function(f)
fpred
fpred(age=30, sex=levels(sex))

anova(f)

p <- Predict(f, age, sex, conf.int=FALSE)
plot(p, ylim=c(-3,3), data=data.frame(age,sex))
# Specifying data to plot.Predict results in sex-specific
# rug plots for age using the Hmisc scat1d function

plsmo(age, L, group=sex, method='raw', add=TRUE, prefix='True', trim=0)
title('Spline Fits with True Regression Functions')

f.bp <- lrm(dz.bp ~ rcs(age,5)*sex + rcs(sys.bp,5))

p <- Predict(f.bp, age, sys.bp, np=75)
for(method in c('contour','persp','image')) {
  bplot(p, method=method)
  #if(method=='image') iLegend(p, c(34,40),c(115, 120))
}


cat('Doing 25 bootstrap repetitions to validate model\n')
validate(f, B=25)   # in practice try to use 300

cat('Doing 25 bootstrap reps to check model calibration\n')
cal <- calibrate(f, B=25)   # use 300 in practice
plot(cal)
title('Calibration of Unpenalized Model')

p <- pentrace(f, penalty=c(.009,.009903,.02,.2,.5,1))

f <- update(f, penalty=p$penalty)
f
specs(f,long=TRUE)
edf <- effective.df(f)

p <- Predict(f, age, sex, conf.int=FALSE)
plot(p, ylim=c(-3,3), data=llist(age, sex))

plsmo(age, L, group=sex, method='raw', add=TRUE, prefix='True', trim=0)
title('Penalized Spline Fits with True Regression Functions')

options(digits=3)
s <- summary(f)
s
plot(s)

s <- summary(f, sex='male')
plot(s)

fpred <- Function(f)
fpred
fpred(age=30, sex=levels(sex))
sascode(fpred)

cat('Doing 40 bootstrap reps to validate penalized model\n')
validate(f, B=40)

cat('Doing 40 bootstrap reps to check penalized model calibration\n')
cal <- calibrate(f, B=40)
plot(cal)
title('Calibration of Penalized Model')

nom <- nomogram(f.bp, fun=plogis,
                funlabel='Prob(dz)',
                fun.at=c(.15,.2,.3,.4,.5,.6,.7,.8,.9,.95,.975))
plot(nom, fun.side=c(1,3,1,3,1,3,1,3,1,3,1))
options(datadist=NULL)

#####################
#Detailed Example 2 #
#####################
# Simulate the data.  
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
treat <- factor(sample(c('a','b','c'), n, TRUE))
num.diseases <- sample(0:4, n, TRUE)
age <- rnorm(n, 50, 10)
cholesterol <- rnorm(n, 200, 25)
weight <- rnorm(n, 150, 20)
sex <- factor(sample(c('female','male'), n, TRUE))
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
cholesterol[1:3] <- NA   # 3 missings, at random

ddist <- datadist(cholesterol, treat, num.diseases,
                  age, weight, sex)
# Could have used ddist <- datadist(data.frame.name)
options(datadist="ddist") # defines data dist. to rms
cholesterol <- impute(cholesterol) # see impute in Hmisc package
# impute, describe, and several other basic functions are
# distributed as part of the Hmisc package


fit <- lrm(y ~ treat*log(cholesterol - 10) +
           scored(num.diseases) +  rcs(age))


describe(y ~ treat + scored(num.diseases) + rcs(age))
# or use describe(formula(fit)) for all variables used in fit
# describe function (in Hmisc) gets simple statistics on variables
#fit <- robcov(fit) # Would make all statistics which follow
                    # use a robust covariance matrix
                    # would need x=TRUE, y=TRUE in lrm
specs(fit) # Describe the design characteristics
a <- anova(fit)
print(a, which='subscripts')          # print which parameters being tested
plot(anova(fit)) # Depict Wald statistics graphically
anova(fit, treat, cholesterol) # Test these 2 by themselves
summary(fit) # Estimate effects using default ranges
plot(summary(fit)) # Graphical display of effects with C.L.
summary(fit, treat="b", age=60) 
# Specify reference cell and adjustment val


summary(fit, age=c(50,70)) # Estimate effect of increasing age from
                           # 50 to 70
summary(fit, age=c(50,60,70)) # Increase age from 50 to 70, 
                              # adjust to 60 when estimating 
                              # effects of other factors
# If had not defined datadist, would have to define
# ranges for all var.


# Estimate and test treatment (b-a) effect averaged
# over 3 cholesterols
contrast(fit, list(treat='b',cholesterol=c(150,200,250)),
              list(treat='a',cholesterol=c(150,200,250)),
         type='average')
# Remove type='average' to get 3 separate contrasts for b-a


# Plot effects.  plot(fit) plots effects of all predictors,
# showing values used for interacting factors as subtitles
# The ref.zero parameter is helpful for showing effects of
# predictors on a common scale for comparison of strength
plot(Predict(fit, ref.zero=TRUE), ylim=c(-2,2))


plot(Predict(fit, age=seq(20,80,length=100), treat, conf.int=FALSE))
# Plots relationship between age and log
# odds, separate curve for each treat, no C.I.
bplot(Predict(fit, age, cholesterol, np=70))
# 3-dimensional perspective plot for age, cholesterol, and
# log odds using default ranges for both variables
p <- Predict(fit, num.diseases, fun=function(x) 1/(1+exp(-x)),
             conf.int=.9)  #or fun=plogis
plot(p, ylab="Prob", conf.int=.9, nlevels=5)
# Treat as categorical variable even though numeric

# Plot estimated probabilities instead of log odds
# Again, if no datadist were defined, would have to
# tell plot all limits
logit <- predict(fit, expand.grid(treat="b",num.diseases=1:3,
                 age=c(20,40,60),
                 cholesterol=seq(100,300,length=10)))
# Also see Predict
#logit <- predict(fit, gendata(fit, nobs=12))
# Interactively specify 12 predictor combinations using UNIX
# For UNIX or Windows, generate 9 combinations with other variables
# set to defaults, get predicted values
logit <- predict(fit, gendata(fit, age=c(20,40,60),
                 treat=c('a','b','c')))


# Since age doesn't interact with anything, we can quickly and
# interactively try various transformations of age,
# taking the spline function of age as the gold standard. We are
# seeking a linearizing transformation.  Here age is linear in the
# population so this is not very productive.  Also, if we simplify the
# model the total degrees of freedom will be too small and
# confidence limits too narrow


ag <- 10:80
logit <- predict(fit, expand.grid(treat="a",
                 num.diseases=0, age=ag,
                 cholesterol=median(cholesterol)),
                 type="terms")[,"age"]
# Also see Predict
# Note: if age interacted with anything, this would be the age
#		"main effect" ignoring interaction terms
# Could also use
#   logit <- plot(f, age=ag, \dots)$x.xbeta[,2]
# which allows evaluation of the shape for any level
# of interacting factors.  When age does not interact with
# anything, the result from
# predict(f, \dots, type="terms") would equal the result from
# plot if all other terms were ignored
# Could also use
#   logit <- predict(fit, gendata(fit, age=ag, cholesterol=median\dots))


plot(ag^.5, logit)  # try square root vs. spline transform.
plot(ag^1.5, logit) # try 1.5 power


# w <- latex(fit)  # invokes latex.lrm, creates fit.tex
# print(w)         # display or print model on screen


# Draw a nomogram for the model fit
plot(nomogram(fit, fun=plogis, funlabel="Prob[Y=1]"))


# Compose S function to evaluate linear predictors from fit
g <- Function(fit)
g(treat='b', cholesterol=260, age=50)
# Leave num.diseases at reference value


# Use the Hmisc dataRep function to summarize sample
# sizes for subjects as cross-classified on 2 key
# predictors
drep <- dataRep(~ roundN(age,10) + num.diseases)
print(drep, long=TRUE)

# Some approaches to making a plot showing how
# predicted values vary with a continuous predictor
# on the x-axis, with two other predictors varying

fit <- lrm(y ~ log(cholesterol - 10) + 
           num.diseases + rcs(age) + rcs(weight) + sex)


combos <- gendata(fit, age=10:100,
                  cholesterol=c(170,200,230),
                  weight=c(150,200,250))
# num.diseases, sex not specified -> set to mode
# can also used expand.grid or Predict


combos$pred <- predict(fit, combos)
require(lattice)
xyplot(pred ~ age | cholesterol*weight, data=combos, type='l')
xYplot(pred ~ age | cholesterol, groups=weight,
       data=combos, type='l') # in Hmisc
xYplot(pred ~ age, groups=interaction(cholesterol,weight),
       data=combos, type='l')


# Can also do this with plot.Predict but a single
# plot may be busy:

ch <- c(170, 200, 230)
p <- Predict(fit, age, cholesterol=ch, weight=150,
     conf.int=FALSE)
plot(p, ~age | cholesterol)

#Here we use plot.Predict to make 9 separate plots, with CLs
p <- Predict(fit, age, cholesterol=c(170,200,230), weight=c(150,200,250))
plot(p, ~age | cholesterol*weight)

options(datadist=NULL)

######################
# Detailed Example 3 #
######################
n <- 2000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
              rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
t <- -log(runif(n))/h
label(t) <- 'Follow-up Time'
e <- ifelse(t<=cens,1,0)
t <- pmin(t, cens)
units(t) <- "Year"
age.dec <- cut2(age, g=10, levels.mean=TRUE)
dd <- datadist(age, sex, age.dec)
options(datadist='dd')
Srv <- Surv(t,e)


# Fit a model that doesn't assume anything except
# that deciles are adequate representations of age
f <- cph(Srv ~ strat(age.dec)+strat(sex), surv=TRUE)
# surv=TRUE speeds up computations, and confidence limits when
# there are no covariables are still accurate.


# Plot log(-log 3-year survival probability) vs. mean age
# within age deciles and vs. sex
p <- Predict(f, age.dec, sex, time=3, loglog=TRUE)
plot(p)
plot(p, ~ as.numeric(as.character(age.dec)) | sex, ylim=c(-5,-1))

# Show confidence bars instead.  Note some limits are not present (infinite)
agen <- as.numeric(as.character(p$age.dec))
xYplot(Cbind(yhat, lower, upper) ~ agen | sex, data=p)



# Fit a model assuming proportional hazards for age and
# absence of age x sex interaction
f <- cph(Srv ~ rcs(age,4)+strat(sex), surv=TRUE)
survplot(f, sex, n.risk=TRUE)
# Add ,age=60 after sex to tell survplot use age=60
# Validate measures of model performance using the bootstrap
# First must add data (design matrix and Srv) to fit object
f <- update(f, x=TRUE, y=TRUE)
validate(f, B=10, dxy=TRUE, u=5)  # use t=5 for Dxy (only)
# Use B=300 in practice
# Validate model for accuracy of predicting survival at t=1
# Get Kaplan-Meier estimates by divided subjects into groups
# of size 200 (for other values of u must put time.inc=u in
# call to cph)
cal <- calibrate(f, B=10, u=1, m=200)  # B=300 in practice
plot(cal)
# Check proportional hazards assumption for age terms
z <- cox.zph(f, 'identity')
print(z); plot(z)


# Re-fit this model without storing underlying survival
# curves for reference groups, but storing raw data with
# the fit (could also use f <- update(f, surv=FALSE, x=TRUE, y=TRUE))
f <- cph(Srv ~ rcs(age,4)+strat(sex), x=TRUE, y=TRUE) 
# Get accurate C.L. for any age
# Note: for evaluating shape of regression, we would not ordinarily
# bother to get 3-year survival probabilities - would just use X * beta
# We do so here to use same scale as nonparametric estimates
f
anova(f)
ages <- seq(20, 80, by=4)   # Evaluate at fewer points. Default is 100
                            # For exact C.L. formula n=100 -> much memory
plot(Predict(f, age=ages, sex, time=3, loglog=TRUE), ylim=c(-5,-1))


# Fit a model assuming proportional hazards for age but
# allowing for general interaction between age and sex
f <- cph(Srv ~ rcs(age,4)*strat(sex), x=TRUE, y=TRUE)
anova(f)
ages <- seq(20, 80, by=6)   
# Still fewer points - more parameters in model


# Plot 3-year survival probability (log-log and untransformed)
# vs. age and sex, obtaining accurate confidence limits
plot(Predict(f, age=ages, sex, time=3, loglog=TRUE), ylim=c(-5,-1))
plot(Predict(f, age=ages, sex, time=3))
# Having x=TRUE, y=TRUE in fit also allows computation of influence stats
r <- resid(f, "dfbetas")
which.influence(f)
# Use survest to estimate 3-year survival probability and
# confidence limits for selected subjects
survest(f, expand.grid(age=c(20,40,60), sex=c('Female','Male')),
        times=c(2,4,6), conf.int=.95)


# Create an S function srv that computes fitted
# survival probabilities on demand, for non-interaction model
f <- cph(Srv ~ rcs(age,4)+strat(sex), surv=TRUE)
srv <- Survival(f)
# Define functions to compute 3-year estimates as a function of
# the linear predictors (X*Beta)
surv.f <- function(lp) srv(3, lp, stratum="sex=Female")
surv.m <- function(lp) srv(3, lp, stratum="sex=Male")
# Create a function that computes quantiles of survival time
# on demand
quant <- Quantile(f)
# Define functions to compute median survival time
med.f <- function(lp) quant(.5, lp, stratum="sex=Female")
med.m <- function(lp) quant(.5, lp, stratum="sex=Male")
# Draw a nomogram to compute several types of predicted values
plot(nomogram(f, fun=list(surv.m, surv.f, med.m, med.f),
         funlabel=c("S(3 | Male)","S(3 | Female)",
                    "Median (Male)","Median (Female)"),
         fun.at=list(c(.8,.9,.95,.98,.99),c(.1,.3,.5,.7,.8,.9,.95,.98),
                   c(8,12),c(1,2,4,8,12))))
options(datadist=NULL)

########################################################
# Simple examples using small datasets for checking    #
# calculations across different systems in which random#
# number generators cannot be synchronized.            #
########################################################

x1 <- 1:20
x2 <- abs(x1-10)
x3 <- factor(rep(0:2,length.out=20))
y  <- c(rep(0:1,8),1,1,1,1)
dd <- datadist(x1,x2,x3)
options(datadist='dd')
f  <- lrm(y ~ rcs(x1,3) + x2 + x3)
f
specs(f, TRUE)
anova(f)
anova(f, x1, x2)
plot(anova(f))
s <- summary(f)
s
plot(s, log=TRUE)
par(mfrow=c(2,2))
plot(Predict(f))
par(mfrow=c(1,1))
plot(nomogram(f))
g <- Function(f)
g(11,7,'1')
contrast(f, list(x1=11,x2=7,x3='1'), list(x1=10,x2=6,x3='2'))
fastbw(f)
gendata(f, x1=1:5)
# w <- latex(f)

f <- update(f, x=TRUE,y=TRUE)
which.influence(f)
residuals(f,'gof')
robcov(f)$var
validate(f, B=10)
cal <- calibrate(f, B=10)
plot(cal)

f <- ols(y ~ rcs(x1,3) + x2 + x3, x=TRUE, y=TRUE)
anova(f)
anova(f, x1, x2)
plot(anova(f))
s <- summary(f)
s
plot(s, log=TRUE)
plot(Predict(f))
plot(nomogram(f))
g <- Function(f)
g(11,7,'1')
contrast(f, list(x1=11,x2=7,x3='1'), list(x1=10,x2=6,x3='2'))
fastbw(f)
gendata(f, x1=1:5)
# w <- latex(f)

f <- update(f, x=TRUE,y=TRUE)
which.influence(f)
residuals(f,'dfbetas')
robcov(f)$var
validate(f, B=10)
cal <- calibrate(f, B=10)
plot(cal)

S <- Surv(c(1,4,2,3,5,8,6,7,20,18,19,9,12,10,11,13,16,14,15,17))
survplot(survfit(S ~ x3))
f <- psm(S ~ rcs(x1,3)+x2+x3, x=TRUE,y=TRUE)
f
# NOTE: LR chi-sq of 39.67 disagrees with that from old survreg
# and old psm (77.65); suspect were also testing sigma=1

for(w in c('survival','hazard'))
 print(survest(f, data.frame(x1=7,x2=3,x3='1'), 
       times=c(5,7), conf.int=.95, what=w))
# S-Plus 2000 using old survival package:
#  S(t):.925 .684 SE:0.729 0.556 Hazard:0.0734 0.255

plot(Predict(f, x1, time=5))
f$var
set.seed(3)
# robcov(f)$var when score residuals implemented
bootcov(f, B=30)$var
validate(f, B=10)
cal <- calibrate(f, cmethod='KM', u=5, B=10, m=10)
plot(cal)
r <- resid(f)
survplot(r)

f <- cph(S ~ rcs(x1,3)+x2+x3, x=TRUE,y=TRUE,surv=TRUE,time.inc=5)
f
plot(Predict(f, x1, time=5))
robcov(f)$var
bootcov(f, B=10)
validate(f, B=10)
cal <- calibrate(f, cmethod='KM', u=5, B=10, m=10)
survplot(f, x1=c(2,19))
options(datadist=NULL)
