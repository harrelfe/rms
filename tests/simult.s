require(rms)
require(multcomp)
set.seed(13)
n <- 200
x1 <- runif(n)
y <- ifelse(runif(n) <= plogis(2*(x1-.5)), 1, 0)
lrm(y ~ x1)
f <- lrm(y ~ rcs(x1,4), x=TRUE, y=TRUE)
g <- bootcov(f, B=1000, coef.reps=TRUE)
anova(f)
specs(f)
# Get simultaneous confidence intervals for estimates at 3 x's
pd <- function(xs) cbind(1, predict(f, data.frame(x1=xs), type='x'))
X <- pd(c(0.05, 0.50, 0.7))
confint(glht(f, X))
# Add a redundant point that does not involve new parameters
X <- pd(c(0.05, 0.50, 0.51, 0.7))
confint(glht(f, X))  # some differences, but slight
# Add a point in a new X space (beyond outer knot)
X <- pd(c(.05, 0.5, 0.51, 0.7, 1))
confint(glht(f, X))

# Add a long sequence of redundant interior points
X <- pd(c(.05, seq(.5, .6, length=100), .7, 1))
confint(glht(f, X))

dd <- datadist(x1); options(datadist='dd')
xs <- seq(0, 1, by=0.02)
i <- Predict(f, x1=xs)
s <- Predict(f, x1=xs, conf.type='simultaneous')
boot <- Predict(g, x1=xs)
b <- rbind(simultaneous=s, individual=i, bootstrap=boot)
plot(b, ~ x1 | .set.)
xYplot(Cbind(yhat,lower,upper) ~ x1, groups=.set., data=b,
       method='bands', type='l', label.curves=list(keys='lines'),
       keyloc=list(x=.1,y=1.5))

contrast(f, list(x1=.2), list(x1=.6))
contrast(f, list(x1=.2), list(x1=c(.6,.8)), conf.type='simult')
