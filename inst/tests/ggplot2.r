## This runs all of the examples in the help file for ggplot2.Predict,
## removing \dontrun{} and not commenting out commands
require(rms)
n <- 500     # define sample size
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

# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male')) +
  .01 * (blood.pressure - 120)
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)
an <- anova(fit)
p <- Predict(fit, age)
ggplot(p)
ggplot(p, ylim=c(-.5, .5))
ggplot(p, xlim=c(40,50))   # problem reported by JoAnn Alvarez

p <- Predict(fit)
# Plot effects in two vertical sub-panels with continuous predictors on top
ggplot(p, sepdiscrete='vertical')
# Plot effects of all 4 predictors with test statistics from anova, and P
ggplot(p, anova=an, pval=TRUE)
ggplot(p, rdata=llist(blood.pressure, age))
# spike histogram plot for two of the predictors
ggplot(p, rdata=llist(blood.pressure, age), sepdiscrete='vertical')

p <- Predict(fit, name=c('age','cholesterol'))   # Make 2 plots
ggplot(p)

p <- Predict(fit, age=seq(20,80,length=100), sex, conf.int=FALSE)
#                        # Plot relationship between age and log
                         # odds, separate curve for each sex,
ggplot(p, subset=sex=='female' | age > 30)
# No confidence interval, suppress estimates for males <= 30

p <- Predict(fit, age, sex)
ggplot(p, rdata=llist(age,sex))
                         # rdata= allows rug plots (1-dimensional scatterplots)
                         # on each sex's curve, with sex-
                         # specific density of age
                         # If data were in data frame could have used that
p <- Predict(fit, age=seq(20,80,length=100), sex='male', fun=plogis)
                         # works if datadist not used
ggplot(p, ylab=expression(hat(P)))
                         # plot predicted probability in place of log odds
per <- function(x, y) x >= 30
ggplot(p, perim=per)       # suppress output for age < 30 but leave scale alone

# Do ggplot2 faceting a few different ways
p <- Predict(fit, age, sex, blood.pressure=c(120,140,160),
             cholesterol=c(180,200,215))
ggplot(p)
ggplot(p, cholesterol ~ blood.pressure)
ggplot(p, ~ cholesterol + blood.pressure)
# color for sex, line type for blood.pressure:
ggplot(p, groups=c('sex', 'blood.pressure'))
# Add legend.position='top' to allow wider plot
# Map blood.pressure to line thickness instead of line type:
ggplot(p, groups=c('sex', 'blood.pressure'), aestype=c('color', 'size'))

# Plot the age effect as an odds ratio
# comparing the age shown on the x-axis to age=30 years

ddist$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: ddist$limits["Adjust to","age"] <- 30
fit <- update(fit)   # make new reference value take effect
p <- Predict(fit, age, ref.zero=TRUE, fun=exp)
# The following fails because of the 3rd element of addlayer
# if ylim is not given, as y=0 is included and can't take log
ggplot(p, ylab='Age=x:Age=30 Odds Ratio', ylim=c(.5, 10),
       addlayer=geom_hline(yintercept=1, col=gray(.7)) +
                geom_vline(xintercept=30, col=gray(.7)) +
                scale_y_continuous(trans='log',
                      breaks=c(.5, 1, 2, 4, 8)))

# Compute predictions for three predictors, with superpositioning or
# conditioning on sex, combined into one graph

p1 <- Predict(fit, age, sex)
p2 <- Predict(fit, cholesterol, sex)
p3 <- Predict(fit, blood.pressure, sex)
p <- rbind(age=p1, cholesterol=p2, blood.pressure=p3)
ggplot(p, groups='sex', varypred=TRUE, adj.subtitle=FALSE, colfill='blue')
ggplot(p, groups='sex', varypred=TRUE, adj.subtitle=FALSE,
       sepdiscrete='vert', rdata=data.frame(age, cholesterol, sex))

# For males at the median blood pressure and cholesterol, plot 3 types
# of confidence intervals for the probability on one plot, for varying age
ages <- seq(20, 80, length=100)
p1 <- Predict(fit, age=ages, sex='male', fun=plogis)  # standard pointwise
p2 <- Predict(fit, age=ages, sex='male', fun=plogis,
              conf.type='simultaneous')               # simultaneous
p3 <- Predict(fit, age=c(60,65,70), sex='male', fun=plogis,
              conf.type='simultaneous')               # simultaneous 3 pts
# The previous only adjusts for a multiplicity of 3 points instead of 100
f <- update(fit, x=TRUE, y=TRUE)
g <- bootcov(f, B=500, coef.reps=TRUE)
p4 <- Predict(g, age=ages, sex='male', fun=plogis)    # bootstrap percentile
p <- rbind(Pointwise=p1, 'Simultaneous 100 ages'=p2,
           'Simultaneous     3 ages'=p3, 'Bootstrap nonparametric'=p4)
# as.data.frame so will call built-in ggplot
ggplot(as.data.frame(p), aes(x=age, y=yhat)) + geom_line() +
 geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2, linetype=0)+
 facet_wrap(~ .set., ncol=2)

# Plots for a parametric survival model
n <- 1000
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
ddist <- datadist(age, sex)
Srv <- Surv(t,e)

# Fit log-normal survival model and plot median survival time vs. age
f <- psm(Srv ~ rcs(age), dist='lognormal')
med <- Quantile(f)       # Creates function to compute quantiles
                         # (median by default)
p <- Predict(f, age, fun=function(x) med(lp=x))
ggplot(p, ylab="Median Survival Time")
# Note: confidence intervals from this method are approximate since
# they don't take into account estimation of scale parameter


# Fit an ols model to log(y) and plot the relationship between x1
# and the predicted mean(y) on the original scale without assuming
# normality of residuals; use the smearing estimator
# See help file for rbind.Predict for a method of showing two
# types of confidence intervals simultaneously.
set.seed(1)
x1 <- runif(300)
x2 <- runif(300)
ddist <- datadist(x1,x2)
y  <- exp(x1+x2-1+rnorm(300))
f <- ols(log(y) ~ pol(x1,2)+x2)
r <- resid(f)
smean <- function(yhat)smearingEst(yhat, exp, res, statistic='mean')
formals(smean) <- list(yhat=numeric(0), res=r[!is.na(r)])
#smean$res <- r[!is.na(r)]   # define default res argument to function
ggplot(Predict(f, x1, fun=smean), ylab='Predicted Mean on y-scale')


# Make an 'interaction plot', forcing the x-axis variable to be
# plotted at integer values but labeled with category levels
n <- 100
set.seed(1)
gender <- c(rep('male', n), rep('female',n))
m <- sample(c('a','b'), 2*n, TRUE)
d <-  datadist(gender, m); options(datadist='d')
anxiety <- runif(2*n) + .2*(gender=='female') + .4*(gender=='female' & m=='b')
tapply(anxiety, llist(gender,m), mean)
f <- ols(anxiety ~ gender*m)
p <- Predict(f, gender, m)
ggplot(p)     # horizontal dot chart; usually preferred for categorical predictors
ggplot(p, flipxdiscrete=FALSE)  # back to vertical
ggplot(p, groups='gender')
ggplot(p, ~ m, groups=FALSE, flipxdiscrete=FALSE)


# Example from Yonghao Pua <puayonghao@gmail.com>
n <- 500  
set.seed(17)  
age <- rnorm(n, 50, 10)
# create an ordinal variable  
duration <- factor(sample(c('None', '10', '20', '30' ,'>30'), n,TRUE))
duration <- factor(duration, levels(duration)[c(5, 1:4)])
# arrange factor levels in ascending order
levels(duration) # shows the correct order "None" "10"   "20"   "30"   ">30"  
label(age) <- 'Age'
label(duration) <- 'Duration'

L <-.045*(age-50) +.01*(duration=='10') +.2*(duration=='20')+
  .3*(duration=='30')+ .9*(duration=='>30')

y <- ifelse(runif(n) < plogis(L), 1, 0)
ddist <- datadist(age, duration)
options(datadist='ddist')
fit <- lrm(y ~ age + duration)
p <- Predict(fit, fun=plogis)
ggplot(p)
ggplot(p, sepdiscrete='vertical', colfill='green', anova=anova(fit))


## From JoAnn Alvarez 2016-10

n <- 800     # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
eyecolor       <- factor(sample(c('green','blue'), n,TRUE))
L <- .4*(sex=='male') +
  .045*(age-50) +
  3*(eyecolor == 'blue')*(sex=='female') +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
y <- ifelse(runif(n) < plogis(L), 1, 0)
ddist <- datadist(age, eyecolor, cholesterol, sex)
options(datadist='ddist')
fit <- lrm(y ~ sex * (eyecolor + age + rcs(cholesterol,4)))
p <- Predict(fit, cholesterol, sex, eyecolor)

ggplot(p)

# Confidence bands automatically suppessed:
ggplot(p, groups = c('eyecolor', 'sex'), aestype=c('color', 'linetype'))

colorscale <- function(...)
  scale_color_manual(...,
                     values=c("#000000", "#E69F00", "#56B4E9",
                              "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
 
ggplot(as.data.frame(p), aes(x=cholesterol, y=yhat, color=eyecolor, linetype=sex)) +
  labs(x=expression(cholesterol), y="log odds",
       title="Adjusted to:age=50.2 ")  +
  geom_line(data=p, mapping=aes(color=eyecolor, linetype=sex)) +
  colorscale(name=expression(eyecolor)) +
  scale_linetype_discrete(name=expression(sex)) +
  theme(legend.position='right') +
#  geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2,
#              linetype=0, fill=I('black'), show.legend=FALSE) +
  coord_cartesian(ylim=c(-4, 6)) +
  theme(plot.title=element_text(size=8, hjust=1))

