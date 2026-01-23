require(rms)
x <- runif(10)
y <- rep(0:4, 2)
f <- orm(y ~ x)
f$coefficients <- c(qlogis(c(.8, .6, .4, .2)), 0)
lp <- predict(f, data.frame(x=0))
lp
plogis(lp)
g <- ExProb(f)
g(lp)
for(i in 0:4) cat('Y>=', i, '', g(lp, y=i), '\n')

set.seed(1)
n     = 200
beta1 = 0.5
beta2 = 0.5
beta3 = 0.3
lambda = 0.5

x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.5)
x3 <- rlnorm(n, 0, 0.6)

ps_true <- plogis(-0.3 + 0.5*x1 + 0.3*x2 + 0.2*x3)
A <- rbinom(n, 1, ps_true)

lin0 <- beta1*x1 + beta2*x2 + beta3*x3
lin1 <- lin0 + lambda

Y <- ifelse(A == 1,
            exp(lin1 + rnorm(n)),
            exp(lin0 + rnorm(n)))
Y <- 1e-7 * round(Y * 1e7)

f <- orm(Y ~ A + x1 + x2 + x3, family = "probit")
g <- ExProb(f)

yu  <- sort(unique(Y))
new <- data.frame(x1, x2, x3, A=1)
lp  <- predict(f, newdata = new)
p1  <- g(lp[1], y=yu)
p   <- g(lp, y = yu)
coef(f)[1:10]
lp[1:10]
p$y[1:10]
p1$prob[1:10]
p$prob[1, 1:10]

