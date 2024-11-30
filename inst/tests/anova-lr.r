require(rms)

## Hauck-Donner effect
set.seed(3)

n  <- 200
x1 <- sample(0:1, n, TRUE)
x2 <- sample(0:1, n, TRUE)
L  <- 0.4 * x1 + 25 * x2
y  <- ifelse(runif(n) <= plogis(L), 1, 0)
# f <- glm.fit(cbind(x1, x2), y, family=binomial())
# f <- glm(y ~ x1 + x2, family=binomial)
# v <- - crossprod(qr.R(f$qr))   #Hessian
# solve(-v, tol=1e-9)
# f <- lrm(y ~ x1 + x2, compvar=FALSE)  # works
f  <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE)
coef(f)

g <- Glm(y ~ x1 + x2, family=binomial, x=TRUE, y=TRUE)
rbind(coef(f), coef(g))
anova(g)
anova(g, test='LR')

x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y  <- round(10*(x1+x2+x3) + 6*rnorm(n))

f <- y ~ x1 + pol(x2, 2) * pol(x3, 2)
s <- Surv(exp(y)) ~ x1 + pol(x2, 2) * pol(x3, 2)
g <- lrm(f, x=TRUE, y=TRUE)
w <- function(fit) {
  print(system.time(a1 <- anova(fit)))
  print(system.time(a2 <- anova(fit, test='LR')))
  print(a1)
  print(a2)
  invisible()
  }
w(g)

g <- orm(f, x=TRUE, y=TRUE)
w(g)

g <- cph(s, x=TRUE, y=TRUE)
w(g)

# survreg.fit2 <- rms:::survreg.fit2
g <- psm(s, x=TRUE, y=TRUE)
w(g)

