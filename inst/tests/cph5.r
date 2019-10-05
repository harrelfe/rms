## Check that the median is NA when there is lots of censoring
require(rms)
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, 
                     rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
dd <- datadist(age, sex)
options(datadist='dd')
S <- Surv(dt,e)

f <- cph(S ~ age + strat(sex), surv=TRUE)
g <- cph(S ~ age + sex, surv=TRUE)
h <- f  # set to f for strat(sex), g for covariate adj for sex
par(mfrow=c(2,2))
for(a in c(45, 55, 65, 75)) {
  survplot(h, sex, age=a)
  title(paste0('Age=', a))
  abline(h=0.5, col=gray(.85))
}

ggplot(Predict(h, age=20:80, sex, time=11))
ggplot(Predict(h, age=20:80, sex, time=12))
ggplot(Predict(h, age=20:80, sex, time=13))
ggplot(Predict(h, age=20:80, sex, time=14))

quan <- Quantile(h)
med  <- function(x) quan(lp=x, stratum=2)
ages <- 70:80
Predict(h, age=ages, sex='Male')$yhat
lp <- predict(h, data.frame(sex='Male', age=ages))
data.frame(age=ages, lp=lp, median=med(lp))

## Estimate survival curves at these medians
if(length(h$strata)) {
  times <- h$time[['sex=Male']]
  surv  <- h$surv[['sex=Male']]
} else {
  times <- h$time
  surv  <- h$surv
}
for(l in lp) print(approx(times, surv ^ exp(l), xout=med(l))$y)

p <- Predict(h, age=ages, fun=med)
plot(p)



