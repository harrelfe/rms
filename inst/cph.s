require(rms)
n <- 2000
set.seed(1)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- 1 + (runif(n)<=.4)
sex <- factor(sample(c('Male','Female'), n, replace=TRUE, prob=c(.6,.4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex==2))
ft <- -log(runif(n))/h
e <- ifelse(ft<=cens,1,0)
print(table(e))
ft <- pmin(ft, cens)
units(ft) <- "Year"
Srv <- Surv(ft, e)
dd <- datadist(age, sex)
options(datadist="dd")
f <- cph(Srv ~ rcs(age,4)+offset(1*(sex=="Male")), eps=1e-9)
g <- coxph(Srv ~ rcs(age,4)+offset(1*(sex=="Male")))
f; g
summary(f)

# Make sure surv summary works
f <- cph(Srv ~ age, surv='summary')
f$surv.summary

# Check relationship between R2 measure and censoring
n <- 2000
set.seed(3)
age <- 50 + 12*rnorm(n)
sex <- factor(sample(c('female','male'), n, TRUE))
cens <- 15*runif(n)
h <- .02 * exp(.1 * (age - 50) + .8 * (sex == 'male'))
t.uncens <- -log(runif(n))/h
e <- ifelse(t.uncens <= cens, 1, 0)
print(table(e))
ft <- pmin(t.uncens, cens)
f <- cph(Surv(ft, e) ~ age + sex, x=TRUE)
f
S <- var(f$x)


cens <- 40*runif(n)
e <- ifelse(t.uncens <= cens, 1, 0)
print(table(e))
ft <- pmin(t.uncens, cens)
g <- cph(Surv(ft, e) ~ age + sex)
g

cens <- 5*runif(n)
e <- ifelse(t.uncens <= cens, 1, 0)
print(table(e))
ft <- pmin(t.uncens, cens)
i <- cph(Surv(ft, e) ~ age + sex)
i

cens <- 2*runif(n)
e <- ifelse(t.uncens <= cens, 1, 0)
print(table(e))
ft <- pmin(t.uncens, cens)
j <- cph(Surv(ft, e) ~ age + sex)
j

# Compute Kent and O'Quigley rho squared W,A tilde
ko <- function(fit, S) {
  cof <- coef(fit)
  rho <- t(cof) %*% S %*% cof
  prn(rho)
  drop(rho / (rho + 1))
}


ko(f, S); ko(g, S); ko(i, S); ko(j, S)

## Compare with OLS R^2
y <- log(h) + rnorm(n)
f <- ols(y ~ age + sex, x=TRUE)
S <- var(cbind(1, f$x))
ko(f, S)

