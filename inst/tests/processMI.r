require(rms)
require(survival)
set.seed(1)
n  <- 400
n1 <- 300; n2 <- 100
data <- data.frame(outcome=c(rnorm(n1, mean = .052, sd = .005),
                             rnorm(n2, mean = .06, sd = .01)),
                   v2=sample(seq(20,80,5),n,T),
                   v3=sample(seq(60,150,1),n,T),
                   v4=c(rnorm(n1, mean = 80, sd = 10),
                        rnorm(n2, mean = 60, sd = 15)),
                   v5=sample(c('M','M','F'),n,T),
                   v6=c(rnorm(n1, mean = 80, sd = 10),
                        rnorm(n2, mean = 120, sd = 30)))

# checking data
head(data)

# setting datadist
dd <- datadist(data);  options(datadist="dd")

# generating missings
m <- function() sample(1:n, 50, FALSE)

for(v in .q(v2, v3, v4, v5, v6)) data[[v]][m()] <- NA

# Imputing
imp <- aregImpute(~ outcome + v2 + v3 + v4 + v5 + v6, data, n.impute=30)

# fitting with validation done separately for each completed dataset
g <- function(fit) list(validate  = validate(fit, B=50),
                        calibrate = calibrate(fit, B=60))

f <- fit.mult.impute(outcome ~ v6 + v2 + rcs(v3) + v5 * rcs(v4),
                     ols, imp, data=data, fun=g, fitargs=list(x=TRUE, y=TRUE))
f
r <- f$funresults
v <- lapply(r, function(x) x$calibrate)
processMI(f, 'validate')
k <- processMI(f, 'calibrate', nind=3)


n <- 350
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n,
                     rep=TRUE, prob=c(.6, .4)))

cens <- 15*runif(n)
h <- .02*exp(.08*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Years"
age[m()] <- NA
sex[m()] <- NA
d <- data.frame(age, sex, dt, e)
dd <- datadist(d)
options(datadist='dd')
rm(age, sex, dt, e)

f <- cph(Surv(dt, e) ~ age + sex, data=d, x=TRUE, y=TRUE, surv=TRUE, time.inc=5)
f$time.inc
a <- aregImpute(~ dt + e + age + sex, data=d, n.impute=10)
g <- function(fit) list(validate =validate(fit,  B=30),
                        calibrate=calibrate(fit, B=30, u=5))
f <- fit.mult.impute(Surv(dt, e) ~ age + sex, cph, a, data=d,
       fun=g, fitargs=list(x=TRUE, y=TRUE, surv=TRUE, time.inc=5))
v <- lapply(f$funresults, function(x) x$validate)
v
processMI(f, 'validate')
k <- lapply(f$funresults, function(x) x$calibrate)
par(mfrow=c(2,3)); for(i in 1:5) plot(k[[i]])
k <- processMI(f, 'calibrate', nind=3)
plot(k)

g <- function(fit) list(calibrate=calibrate(fit, B=30, u=5, cmethod='KM', m=50))
f <- fit.mult.impute(Surv(dt, e) ~ age + sex, cph, a, data=d,
                     fun=g, fitargs=list(x=TRUE, y=TRUE, surv=TRUE, time.inc=5))
k <- processMI(f, 'calibrate', nind=3)
