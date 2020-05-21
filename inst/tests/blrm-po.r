doit <- FALSE
if(doit) {

require(rms)
stanSet()
set.seed(1)
n  <- 100
x  <- rnorm(n)
x2 <- rnorm(n)
x  <- x  - mean(x)
x2 <- x2 - mean(x2)
Y <- x + 2 * x2 + rnorm(n)

y <- round(3*Y)
i <- rep(1 : 20, 5)
f <- lrm(y ~ x + x2)
b <- blrm(y ~ x + x2, priorsd=1000)
bc <- blrm(y ~ x + x2 + cluster(i), priorsd=1000)
cbind(lrm=coef(f), blrm=coef(b, 'mean'), 'cluster blrm'=coef(bc, 'mean'))

d <- blrm(y ~ x + x2, priorsd=1000, standata=TRUE)
saveRDS(d, '~/tmp/d.rds')

set.seed(1)
n  <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 3 * x1 + 4 * x2 + 2 * rnorm(n)
lrm(cut2(y, g=20) ~ pol(x1,2) + pol(x2,2))

do <- function(fun, eq=TRUE) {
  par(mfrow=c(5,6), mar=c(3, 2.5, 1, 1))
  m <- 0
  gs <- c(2:4, seq(5, 100, by=5))
  for(g in gs) {
    u <- if(eq) cut2(y, g=g)
         else cut2(y, if(g==2) 0 else  seq(-10, 10, length=g-1))
    if(length(unique(u)) != g) next
    f <- lrm (u ~ pol(x1,2) + x2)
    b <- blrm(u ~ pol(x1,2) + x2, conc=fun(g), priorsd=10000)
    u <- coef(f)[1 : (g-1)]
    v <- coef(b)[1 : (g-1)]
    plot(u, v)
    abline(a=0, b=1, col='red')
    d <- abs(u - v)
    w <- paste(g, round(max(d), 3), round(mean(d), 3))
    title(w)
    m <- m + mean(d)
  }
  cat('Avg mean abs diff:', m / length(gs), '\n')
}

do(function(g) 2/max(1, g-5))
do(function(g) 1/g)            # 0.042
do(function(g) 1/(g + 2))
do(function(g) 1/sqrt(g))
do(function(g) 1/(max(3, g)))
do(function(g) 1/min(g, 20))   # not bad esp late 0.026
do(function(g) 1/min(g, 30))   # not as good
do(function(g) 1 / (2 + (g/3)))  # winner below
do(function(g) 1/(0.8 + 0.35 * max(g, 3)), eq=FALSE)  # good w/n=1000
do(function(g) 1/(0.8 + 0.35 * max(g, 3)), eq=TRUE)   # good w/n=1000


set.seed(1)
n  <- 200
x1  <- rnorm(n)
x2 <- rnorm(n)
y <- 3 * x1 + 4 * x2 + 2 * rnorm(n)

z <- expand.grid(g=c(2:5, seq(10, 90, by=5)),
                 conc=c(0.01, 0.02, .03, .04, seq(0.05, 1, by=0.05)))
z <- data.frame(z, mad=NA)
for(i in 1 : nrow(z)) {
  g    <- z$g[i]
  conc <- z$conc[i]
  cat('conc=', conc, '\n')
  u <- cut2(y, g=g)
  f <- lrm(u ~ x1 + x2)
  b <- blrm(u ~ x1 + x2, conc=conc, priorsd=10000)
  u <- coef(f)[1 : (g-1)]
  v <- coef(b)[1 : (g-1)]
  d <- abs(u - v)
  z$mad[i] <- mean(d)
}

saveRDS(z, 'blrm-po.rds') 

ggplot(z, aes(x=g, y=conc, size=log(mad))) + geom_point()

gs <- unique(z$g)
bestconc <- numeric(length(gs))
i <- 0
for(k in gs) {
  i <- i + 1
  u <- subset(z, g==k)
  bestconc[i] <- with(u, conc[which.min(mad)])
}
plot(gs, bestconc)
plot(gs, 1/bestconc)
abline(lsfit(gs, 1/bestconc))
summary(lm(1/bestconc ~ gs))
summary(lm(1/bestconc ~ pol(gs, 2)))

## 1 / (2 + k/3)

1/bestconc[gs==2]
1/bestconc[gs==3]

plot(gs, 1/bestconc)
abline(lsfit(gs[gs > 2], 1/bestconc[gs > 2]))
lm(1/bestconc ~ gs, subset=gs > 2)
}
