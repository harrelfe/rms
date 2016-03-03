# Thanks: Kai Chen. M.D.
# Resident
# Dept of Breast surgery,
# Breast Tumor Center;  
# Sun Yat-sen Memorial Hospital
# chenkai23@mail.sysu.edu.cn

# Fit logistic model on 100 observations simulated from the actual 
# model given by Prob(Y=1 given X1, X2, X3) = 1/(1+exp[-(-1 + 2X1)]),
# where X1 is a random uniform [0,1] variable.  Hence X2 and X3 are 
# irrelevant.  After fitting a linear additive model in X1, X2,
# and X3, the coefficients are used to predict Prob(Y=1) on a
# separate sample of 100 observations.  Note that data splitting is
# an inefficient validation method unless n > 20,000.

require(rms)
set.seed(1)
n <- 200
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
logit <- 2*(x1-.5)
P <- 1/(1+exp(-logit))
y <- ifelse(runif(n)<=P, 1, 0)
d <- data.frame(x1,x2,x3,y)
f <- lrm(y ~ x1 + x2 + x3, subset=1:100)
pred.logit <- predict(f, d[101:200,])
phat <- 1/(1+exp(-pred.logit))

val.prob(phat, y[101:200], m=20, cex=.5)  # subgroups of 20 obs.

val.prob(x2[101:200], y[101:200], m=20, cex=.5)  # subgroups of 20 obs.
