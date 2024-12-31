## Test transx option for lrm
require(rms)
set.seed(3)
x <- rnorm(30)
y <- sample(0:4, 30, TRUE)
h <- orm(y ~ pol(x, 2))
f <- lrm(y ~ pol(x, 2), transx=FALSE)
range(vcov(f, intercepts='all') - vcov(h, intercepts='all')) # correct
g <- lrm(y ~ pol(x, 2), transx=TRUE)
range(coef(f) - coef(g))
range(vcov(f) - vcov(g))
vcov(f) / vcov(g)

f <- orm(y ~ pol(x, 2), scale=FALSE)
g <- orm(y ~ pol(x, 2), scale=TRUE)
fi <- f$info.matrix
gi <- g$info.matrix
# For f info matrix is for original x
# For g info matrix is for centered/scaled x and reverse scaling is
# done by vcov.orm/vcov.lrm using infoMxop
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

