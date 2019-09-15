# Test execution time of cluster sandwich covariance estimator with logistic models
# GEE using a working independence binary logistic regression model
# 1,000,000 records on 100,000 subjects, 10 covariates
# Times are on a Lenovo X1 laptop running Linux
require(rms)
set.seed(1)
n <- 1000000
subjects <- 100000
y <- sample(0:1, n, TRUE)
x <- matrix(runif(10*n), ncol=10)
id <- sample(1:subjects, n, TRUE)
# Fit binary logistic model with working independence structure
system.time(f <- lrm(y ~ x, x=TRUE, y=TRUE))   # 4s
# g will have betas from f but with robust covariance matrix
system.time(g <- robcov(f, id))                # 1.4s
diag(vcov(f)) / diag(vcov(g))

# Check robcov's ability to ignore duplicate data
m <- n / subjects
set.seed(1)
y <- rep(sample(0 : 1, subjects, replace=TRUE), each=m)
id <- rep(1 : subjects, each=m)
j <- 1
x <- matrix(NA, nrow=n, ncol=10)
for(i in 1 : subjects) {
    x[j : (j + m - 1), ] <- matrix(rep(runif(10), each=m), nrow=m)
    j <- j + m
}
f <- lrm(y ~ x, x=TRUE, y=TRUE)
g <- robcov(f, id)
diag(vcov(f) / vcov(g))
