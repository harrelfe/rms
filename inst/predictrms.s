## David van Klaveren
## Erasmus MC, Department of Public Health, room Ae-110
## E-mail   d.vanklaveren.1@erasmusmc.nl
## se.fit comparisons before predictrms changed to use original covariate
## means instead of "adjust to" values for cph

require(rms)
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c("Male","Female"), n, rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=="Female"))
dt <- -log(runif(n))/h
label(dt) <- "Follow-up Time"
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
dd <- datadist(age, sex)
options(datadist="dd")
Srv <- Surv(dt,e)

f <- cph(Srv ~ sex + age , x=TRUE, y=TRUE, se.fit=TRUE)
## skipped splines to keep example simple
p <- predict(f,newdata=data.frame(sex,age), se.fit=TRUE)
## predict with newdata = original data

## linear predictors are equal for cph and predict:
sum((f$linear.predictors-p$linear.predictors)^2)

## and so are se after predictrms fixed
sum((f$se.fit-p$se.fit)^2)

### Reconstruction of difference
X    <- f$x
beta <- f$coef
cov  <- f$var

X.center.mean   <- sweep(X, 2, c(mean(sex=="Male"),  mean(age)))
X.center.median <- sweep(X, 2, c(median(sex=="Male"),median(age)))

lp.center.mean <- X.center.mean%*%beta
se.center.mean <- drop(sqrt(((X.center.mean %*% cov) * X.center.mean) %*%
                            rep(1, ncol(X.center.mean))))
se.center.median <- drop(sqrt(((X.center.median %*% cov) * X.center.median) %*%
                              rep(1, ncol(X.center.median))))

## linear predictors are equal for fit/predict and mean centered lp:
sum((f$linear.predictors-lp.center.mean)^2)

## cph$se.fit is equal to mean centered se
sum((f$se.fit-se.center.mean)^2)

## predict$.se.fit is no longer equal to median centered se
sum((p$se.fit-se.center.median)^2)


## Check ref.zero=TRUE
set.seed(1)
n <- 30
x1 <- 100 + rnorm(n)
x2 <-   5 + rnorm(n)
x3 <- c(rep('a', 5), rep('b', 5), rep('c', 20))
dd <- datadist(x1, x2, x3); options(datadist='dd'); dd
resid <- rnorm(n)
y <- x1 + x2 + .5*(x3 == 'b') + 1*(x3 == 'c') + resid

f <- ols(y ~ pol(x1, 2) + x2 + x3)
f
w <- data.frame(x1=100.4, x2=4.5, x3='c')
predict(f, w, conf.int=.95) 
predict(f, w, type='adjto')
c(median(x1), median(x1)^2, median(x2))
k <- coef(f)
ycenter <- k[1] + k[2]*median(x1) + k[3]*median(x1)^2 + k[4]*median(x2) + k[6]
ycenter
predict(f, w, conf.int=.95, ref.zero=TRUE)
k[1] + k[2]*100.4 + k[3]*100.4^2 + k[4]*4.5 + k[6] - ycenter
