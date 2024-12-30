## Test transx option for lrm
require(rms)
x <- rnorm(30)
y <- sample(0:4, 30, TRUE)
f <- lrm(y ~ pol(x, 2), transx=FALSE)
g <- lrm(y ~ pol(x, 2), transx=TRUE)
coef(f) - coef(g)
d <- vcov(f) - vcov(g)  # way off
d
max(abs(d))

f <- orm(y ~ pol(x, 2), scale=FALSE)
g <- orm(y ~ pol(x, 2), scale=TRUE)
fi <- f$info.matrix
gi <- g$info.matrix
for(n in names(fi)) {
  cat(n, '\n')
  print(fi[[n]])
  print(gi[[n]])
}

coef(f)
coef(g)
vcov(f) / vcov(g)
vcov(f, regcoef.only=FALSE) / vcov(g, regcoef.only=FALSE)
vcov(f, intercepts='all') / vcov(g, intercepts='all')
vcov(f, intercepts=1) / vcov(g, intercepts=1)
vcov(f, intercepts=c(1,3)) / vcov(g, intercepts=c(1,3))
