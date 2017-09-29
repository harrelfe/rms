require(rms)
set.seed(1)
n <- 20
x <- as.matrix(1:n)
#x <- cbind(1:n, (1:n)^2)
#colnames(x) <- 'age'
y <- sample(0:1, n, TRUE)
f <- lrm(y ~ x)


N <- 100
set.seed(1)
time <- rexp(N)
status <- sample(0:1, N, replace = TRUE)
S <- Surv(time, status)
x1 <- gl(2, 50)
x2 <- runif(N)
x3 <- sample(1:3, N, replace=TRUE)

ols(time ~ x1)
ols(time ~ scored(x3))
ols(time ~ catg(x3))


# Makes last colname x1 %ia% x2 which is really inconsistent:
model.matrix(~ x1 + rcs(x2) + x1 %ia% x2)
x3 <- c(rep('A', 33), rep('B', 33), rep('C', 34))
x4 <- runif(N) > 0.5
# Makes last 2 colnames x3 %ia% x2x3=B * x2,  x3 %ia% x2x3=C * x2
model.matrix(~ x3 + rcs(x2) + x3 %ia% x2)
cph(S ~ x3 + rcs(x2) + x3 %ia% x2)

ols(time ~ x1 + rcs(x2) + x1 %ia% x2)
lrm(status ~ x1 + rcs(x2) + x1 %ia% x2)
options(debug=TRUE,width=110)
cph(S ~ x1 + rcs(x2) + x1 %ia% rcs(x2))
cph(S ~ x1 + rcs(x2) + x1 %ia% x2)

cph(S ~ x1 * rcs(x2))

ols(time ~ x1 + x4)
cph(S ~ x1 + x4)
colnames(model.matrix(~ x1 + x4 + x1 %ia% x4))
cph(S ~ x1 + x4 + x1 %ia% x4)


## From https://github.com/harrelfe/rms/issues/29#issuecomment-303423887
## https://github.com/harrelfe/rms/issues/29#issuecomment-328581864
d <- expand.grid(
X1 = factor(c('05: X1 <= 178','01: X1 <= 6', '03: X1 <= 52', '05: X1 <= 178')),
X2 = factor(c('04: X2 <= 75','01: X2 <= 6', '05: X2 > 75', '05: X2 > 75')),
X3 = factor(c('04: X3 <= 552','01: X3 <= 1', '04: X3 <= 552', '06: X3 > 1313')),
rep = 1 : 100)
set.seed(1)
d$TARGET <- sample(0 : 1, nrow(d), replace=TRUE)

lrm(TARGET ~ ., data = d)
options(debug=TRUE)
cph(Surv(TARGET) ~ ., data=d)
