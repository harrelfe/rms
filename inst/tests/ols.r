require(rms)
n = 800
x = rnorm(n)
x[4:6] = NA
y = x + rnorm(n)
fit = ols(y ~ x)
print(fit, latex=TRUE)
