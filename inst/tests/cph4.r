require(rms)
options(debug=TRUE)
set.seed(1)
n  <- 100
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- sample(0:4, n, TRUE)
x5 <- ordered(sample(1:3, n, TRUE), levels=1:3, labels=c('I','II','III'))
S  <- Surv(runif(n))
# f <- cph(S ~ x1 + pol(x2, 2) + rcs(x3, 4) + scored(x4) + x5)   # FAILS
options(contrasts=c('contr.treatment', 'contr.treatment'))
f <- cph(S ~ x1 + pol(x2, 2) + rcs(x3, 4) + scored(x4) + x5)
f
     
