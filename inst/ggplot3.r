# Test anova= in ggplot
require(rms)
set.seed(1)
x1 <- runif(100)
x2 <- runif(100)
x3 <- sample(c('a','b'), 100, TRUE)
x4 <- sample(c('k','l','m'), 100, TRUE)
y <- runif(100)
dd <- datadist(x1, x2, x3, x4); options(datadist='dd')
f <- ols(y ~ x1 + x2 + x3 + x4)

a <- anova(f)
ggplot(Predict(f), anova=a)   # ok
ggplot(Predict(f), anova=a, sepdiscrete='vertical')
