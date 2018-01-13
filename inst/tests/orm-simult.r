require(rms)
set.seed(1)

xa <- runif(100, 0, 5)
xb <- runif(100, 0, 5)

ya <- xa + rnorm(100, 0, 1)
yb <- xb + rnorm(100, 0.25, 1)

arma <- rep("Placebo", 100)
armb <- rep("Intervention", 100)

x <- c(xa, xb)
y <- c(ya, yb)
arm <- c(arma, armb)

dd <- datadist(x, y, arm)
options(datadist="dd")

f_ols <- ols(y ~ x + arm)
#Works
c_ols <- contrast(f_ols, list(arm="Placebo"), list(arm="Intervention"), conf.type='simultaneous')

f_orm <- orm(y ~ x + arm)
#Fails with dimensions of coefficients and covariance matrix don't match 
c_orm <- contrast(f_orm, list(arm="Placebo"), list(arm="Intervention"), conf.type='simultaneous')
#Error in modelparm.default(model, ...) : 
#  dimensions of coefficients and covariance matrix don't match

