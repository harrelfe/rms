require(rms)
tt <- function() {
  dp <- function(...) {
    pdf('/tmp/z.pdf')
    print(...)
    dev.off()
  }
n <- 350     # define sample size
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

ddist <<- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)
an <- anova(fit)
# Plot effects in two vertical sub-panels with continuous predictors on top
# ggplot(Predict(fit), sepdiscrete='vertical')
# Plot effects of all 4 predictors with test statistics from anova, and P
dp(ggplot(Predict(fit), anova=an, pval=TRUE))
# ggplot(Predict(fit), rdata=llist(blood.pressure, age))
# spike histogram plot for two of the predictors

# p <- Predict(fit, name=c('age','cholesterol'))   # Make 2 plots
# ggplot(p)

# p <- Predict(fit, age=seq(20,80,length=100), sex, conf.int=FALSE)
#                        # Plot relationship between age and log
                         # odds, separate curve for each sex,
# ggplot(p, subset=sex=='female' | age > 30)
# No confidence interval, suppress estimates for males <= 30

# p <- Predict(fit, age, sex)
# ggplot(p, rdata=llist(age,sex))
                         # rdata= allows rug plots (1-dimensional scatterplots)
                         # on each sex's curve, with sex-
                         # specific density of age
                         # If data were in data frame could have used that
# p <- Predict(fit, age=seq(20,80,length=100), sex='male', fun=plogis)
                         # works if datadist not used
# ggplot(p, ylab=expression(hat(P)))
                         # plot predicted probability in place of log odds
# per <- function(x, y) x >= 30
# ggplot(p, perim=per)       # suppress output for age < 30 but leave scale alone

# Do ggplot2 faceting a few different ways
p <- Predict(fit, age, sex, blood.pressure=c(120,140,160),
             cholesterol=c(180,200,215))
# ggplot(p)
dp(ggplot(p, cholesterol ~ blood.pressure))
# ggplot(p, ~ cholesterol + blood.pressure)
# color for sex, line type for blood.pressure:
dp(ggplot(p, groups=c('sex', 'blood.pressure')))
# Add legend.position='top' to allow wider plot
# Map blood.pressure to line thickness instead of line type:
# ggplot(p, groups=c('sex', 'blood.pressure'), aestype=c('color', 'size'))

# Plot the age effect as an odds ratio
# comparing the age shown on the x-axis to age=30 years

# ddist$limits$age[2] <- 30    # make 30 the reference value for age
# Could also do: ddist$limits["Adjust to","age"] <- 30
# fit <- update(fit)   # make new reference value take effect
# p <- Predict(fit, age, ref.zero=TRUE, fun=exp)
# ggplot(p, ylab='Age=x:Age=30 Odds Ratio',
#        addlayer=geom_hline(yintercept=1, col=gray(.8)) +
#                 geom_vline(xintercept=30, col=gray(.8)) +
#                 scale_y_continuous(trans='log',
#                       breaks=c(.5, 1, 2, 4, 8))))

# Compute predictions for three predictors, with superpositioning or
# conditioning on sex, combined into one graph

p1 <- Predict(fit, age, sex)
p2 <- Predict(fit, cholesterol, sex)
p3 <- Predict(fit, blood.pressure, sex)
p <- rbind(age=p1, cholesterol=p2, blood.pressure=p3)
dp(ggplot(p, groups='sex', varypred=TRUE, adj.subtitle=FALSE))
# ggplot(p, groups='sex', varypred=TRUE, adj.subtitle=FALSE, sepdiscrete='vert')


# Make an 'interaction plot', forcing the x-axis variable to be
# plotted at integer values but labeled with category levels
n <- 100
set.seed(1)
gender <- c(rep('male', n), rep('female',n))
m <- sample(c('a','b'), 2*n, TRUE)
d <<-  datadist(gender, m); options(datadist='d')
anxiety <- runif(2*n) + .2*(gender=='female') + .4*(gender=='female' & m=='b')
tapply(anxiety, llist(gender,m), mean)
f <- ols(anxiety ~ gender*m)
p <- Predict(f, gender, m)
# ggplot(p)     # horizontal dot chart; usually preferred for categorical predictors
# ggplot(p, flipxdiscrete=FALSE)  # back to vertical
dp(ggplot(p, groups='gender'))
dp(ggplot(p, ~ m, groups=FALSE, flipxdiscrete=FALSE))
}
print(system.time(tt()))

