require(rms)
n <- 1e6
set.seed(8)
y <- sample(1:5, n, prob=c(.1, .2, .35, .35, .05), replace=TRUE)
table(y)
x <- sample(0:1, n, replace=TRUE)
means <- tapply(y, x, mean)
means
dd <- datadist(x); options(datadist='dd')
f <- orm(y ~ x)
M <- Mean(f)
lp <- Predict(f)   # same as Predict(f, x=0:1) here
lp
lp <- lp$yhat
lp
M(lp)
X <- predict(f, data.frame(x=0:1), type='x')
M(lp, X, conf.int=0.95)
