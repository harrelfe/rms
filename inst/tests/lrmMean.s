require(rms)
set.seed(6) # was 3
n <- 85+15
y <- sample(1:10, n, TRUE)
x1 <- runif(n)
f <- lrm(y ~ x1, x=TRUE, y=TRUE)
g <- bootcov(f, B=500, coef.reps=TRUE)
m <- Mean(f)
Predict(f, x1=c(.25, .75), fun=m)
Predict(g, x1=c(.25, .75), fun='mean')
h <- ols(y ~ x1)
Predict(h, x1=c(.25, .75))
