# https://github.com/harrelfe/rms/issues/55

require(rms)
set.seed(1)
n <- 20
X2 <- factor(c(rep(0, n/2), rep(1, n/2)))
X21 <- rep(1 : (n/2), 2)
y <- rep(0 : 1, n/2)
options(rmsdebug=TRUE)
f <- lrm(y ~ X2 + X21, method='model.frame')
attributes(f)$Design$mmcolnames

# Problem is inherent to R
colnames(model.matrix(~ X2 + X21))

## https://github.com/harrelfe/rms/issues/29#issuecomment-417901353
d <- data.frame(
  X = sample(1:700),
  Y = sample(c("yes", "no"),700, replace = TRUE),
  Z = sample (c("Back pain", "Leg Pain", "Back pain = Leg pain"),700, replace = TRUE)
)
options(rmsdebug=TRUE)
lrm(Y~X+Z, data=d)

d$Z <- sample(c("Back pain", "Leg Pain", "Back pain = Leg pain"),700, replace = TRUE)
lrm(Y~X+Z, data=d)

