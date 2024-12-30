require(rms)
n <- 250
set.seed(8)
y <- sample(1:5, n, prob=c(.1, .2, .35, .35, .05), replace=TRUE)
table(y)
x <- sample(0:1, n, replace=TRUE)
wt <- sample(1:3, n, TRUE)
f <- lrm(y ~ x)
coef(f)

f <- lrm(y ~ x, weights=wt)
coef(f)
g <- orm(y ~ x, weights=wt)
coef(g)
coef(f) - coef(g)
vcov(f) / vcov(g, intercepts='all')

infoMxop(g$info.matrix) - infoMxop(f$info.matrix)
infoMxop(g$info.matrix, invert=TRUE) - vcov(f)
