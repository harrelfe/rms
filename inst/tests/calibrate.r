library(rms)
n <- 1000 # define sample size
set.seed(17) # so can reproduce the results
age <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol <- rnorm(n, 200, 25)
sex <- factor(sample(c('female','male'), n,TRUE))
label(age) <- 'Age' # label is in Hmisc
label(cholesterol) <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex) <- 'Sex'
units(cholesterol) <- 'mg/dl' # uses units.default in Hmisc
units(blood.pressure) <- 'mmHg'
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
y <- ifelse(runif(n) < plogis(L), 1, 0)
cholesterol[1:3] <- NA # 3 missings, at random
ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')
f <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
         x=TRUE, y=TRUE)                                                      

cal <- calibrate(f, B=80)
class(cal)

plot(cal, 
     xlim=c(0,1.0),
     ylim=c(0,1.0),
     xlab="Predicted Probability of 1y Event-Free", 
     ylab="Actual 1y",
     subtitles=TRUE,
     riskdist=FALSE, 
     scat1d.opts=list(nhistSpike=200))
