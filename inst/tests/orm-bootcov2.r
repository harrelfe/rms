# From doc/rms/rmsc/validate.qmd

require(rms)
n <- 60
set.seed(3)
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rnorm(n)
x5 <- rnorm(n)
y  <- round(x1 + x2 + rnorm(n), 2)  # used ,1 for n=1200
d  <- data.frame(x1, x2, x3, x4, y)
f  <- orm(y ~ pol(x1, 2) * pol(x2, 2) + x3 + x4 + x5,
          x=TRUE, y=TRUE)
i <- (num.intercepts(f) + 1) : length(coef(f))
iref <- f$interceptRef
jj <- c(iref, i)
Ko <- coef(f)[jj]

g <- lrm(y ~ pol(x1, 2) * pol(x2, 2) + x3 + x4 + x5,
         x=TRUE, y=TRUE)
range(vcov(f, intercepts='all') / vcov(g))

h <- MASS::polr(factor(y) ~ pol(x1, 2) * pol(x2, 2) + x3 + x4 + x5)
coef(f)[i] - coef(h)
v <- vcov(h)

k <- num.intercepts(f)
p <- length(coef(f)) - k
ia <- c((p + 1):(k + p), 1 : p)
ir <- c(p + iref, 1 : p)

diag(vcov(f)) / diag(v[ir, ir])
range(diag(vcov(f, intercepts='all')) / diag(v[ia, ia]))

# Manual bootcov
set.seed(1)
B <- 400
co <- matrix(NA, B, p + 1)
ytarget <- median(y)
ytarget
X <- f$x

a <- 0
for(i in 1 : B) {
  j <- sample(n, n, TRUE)
  b <- orm.fit(X[j, ], y[j])
  if(b$fail) next
  a <- a + 1
  # Keep only the intercept closest to the target
  yu <- b$yunique[-1]
  m  <- which.min(abs(yu - ytarget))
  nint <- num.intercepts(b)
  # cat('i=', i, ' nint=', nint, ' m=', m, ' yu=', yu[m], '\n')
  cof <- coef(b)
  co[a, ] <- cof[c(m, (nint + 1) : length(cof))]
}

co <- co[1 : a, ]
dim(co)
rn <- function(x) round(x, 3)
rn(co[1:20,])
rn(cbind(apply(co, 2, mean), apply(co, 2, median), Ko))

rn(diag(var(co)) / diag(vcov(f)))
# Bootstrap variances 4x larger than MLE estimates with n=60
# with n=1200 bootstrap about 1.3x larger

set.seed(3)
b1 <- bootcov(f, B=400)
vb <- b1$var[jj, jj]
rn(diag(vb) / diag(vcov(f)))   # 4x increase for bootstrap

set.seed(3)
b2 <- bootcov(f, B=400, ytarget=NA)   # save only one intercept
rn(diag(b2$var) / diag(vcov(f)))      # 4x increase for bootstrap

rn(diag(b1$var[jj, jj]) / diag(b2$var))    # ratios all 1.0

rx <- rexVar(b1, data=d)
rx
rx <- rexVar(b2, data=d)
rx
plot(rx)
