require(rms)
set.seed(1)
n <- 20
y <- sample(1:4, n, replace=TRUE)
x1 <- runif(n)
x2 <- runif(n)
d <- data.frame(x1, x2)
f <- lrm(y ~ x1 + x2)
s <- 1:4
f$linear.predictors[s]
options(digits=3)
predict(f)[s]    # kint=1
predict(f, d)[s] # kint=1
xb <- as.vector(cbind(x1, x2) %*% as.matrix(coef(f)[c('x1','x2')]))
coef(f)[1] + xb[s] # kint=1

predict(f, type='fitted.ind')[s,]
# get Prob(Y=4)
plogis(coef(f)[3] + xb)[s]

g <- lrm(y ~ x1 + x2, x=TRUE, linear.predictors=FALSE)
predict(g)[s]  # correct; kint=1
predict(g, d)[s] # agrees (kint=1)
predict(g, type='fitted.ind')[s,]
