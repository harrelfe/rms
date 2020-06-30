## Test partial proportional odds model where the tau matrix is not square
require(rms)
require(VGAM)
stanSet()
set.seed(1)
n <- 300
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
y  <- sample(1:4, n, TRUE)

f <- vgam(y ~ x1 + x2 + x3, cumulative(reverse=TRUE, parallel=FALSE))
coef(f)

# Reparameterize since vgam uses total effects and blrm uses
# increments in log odds due to partial PO effects
b <- blrm(y ~ x1 + x2 + x3, ~ x1 + x2 + x3)
u <- coef(f)[c(1, 2, 3, 4, 7, 10, 5, 8, 11, 6, 9, 12)]
v <- coef(b, 'mode')
v[7] <- v[4] + v[7]
v[8] <- v[5] + v[8]
v[9] <- v[6] + v[9]
v[10] <- v[4] + v[10]
v[11] <- v[5] + v[11]
v[12] <- v[6] + v[11]

cbind(u, v)
