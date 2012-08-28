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
