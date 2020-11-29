require(rms); require(MASS)

set.seed(1)
n <- 100
y <- sample(1:8, n, TRUE)
#y <- runif(n)
x1 <- sample(c(-1,0,1), n, TRUE)
x2 <- sample(c(-1,0,1), n, TRUE)
system.time(f <- lrm(y ~ x1 + x2, eps=1e-5, trace=TRUE))
fkeep <- f

# xless(solve(vcov(f)))
predict(f, data.frame(x1=0,x2=0), type='fitted.ind')

h <- polr(as.factor(y) ~ x1 + x2, Hess=TRUE)
v <- solve(vcov(h))
s <- c(3:ncol(v),1:2)
v[s,s]
v[s,s] / solve(vcov(f))


g <- orm(y ~ x1 + x2, eps=1e-5, trace=TRUE)

system.time(g <- orm(y ~ x1 + x2, eps=.001)) #, trace=TRUE)
coef(g) - coef(f)
w <- vcov(g, intercepts='all') / vcov(f) - 1
max(abs(w))


m <- Mean(g)
formals(m) <- list(lp=NULL, intercepts=runif(30000), values=runif(30001),
                   conf.int=0, interceptRef=3, cumprob=function(x) 1 / (1 + exp(-x)))
system.time(m(1))
system.time(m(1:100))
system.time(m(1:1000))

set.seed(1)
n <- 1000
x1 <- c(rep(0,500), rep(1,500))
y <- rnorm(n) + 3*x1
g <- orm(y ~ x1)
g
k <- coef(g)
plot(1:999, k[1:999])
h <- orm(y ~ x1, family='probit')
plot(coef(g)[1:999], coef(h)[1:999])

tapply(y, x1, mean)

m <- Mean(g)
m(w <- k[g$interceptRef] + k['x1']*c(0,1))
#mf <- Mean(f)
#k <- coef(f)
#mf(k[1] + k['x1']*c(0,1))
mh <- Mean(h)
wh <- coef(h)[h$interceptRef] + coef(h)['x1']*c(0,1)
mh(wh)

qu <- Quantile(g)
qu(.1, w)
qu(.5, w)
qu(.9, w)
tapply(y, x1, quantile, probs=c(.1,.5,.9))

set.seed(1)
n <- 1000
x1 <- c(rep(0,500), rep(1,500))
y <- exp(rnorm(n) + 3*x1)
g <- orm(y ~ x1)
g
k <- coef(g)
plot(1:999, k[1:999])
m <- Mean(g)
m(w <- k[1] + k['x1']*c(0,1))
m(w <- k[g$interceptRef] + k['x1']*c(0,1))
tapply(y, x1, mean)

qu <- Quantile(g)
qu(.1, w)
tapply(y, x1, quantile, probs=.1)
qu(.5, w)
tapply(y, x1, quantile, probs=.5)
qu(.9, w)
tapply(y, x1, quantile, probs=.9)

## Check quantile calculations
qu <- Quantile(g)
## .9 = Prob(Y >= 2) .8 = Prob(Y >= 3) etc.
## Prob Y <= j, j = 1, ... 10 = .1, .2, ..., 1
## .1 quantile = 1, .2 quantile = 2, ..., .9 quantile = 9
formals(qu) <- list(q=.5, lp=0, intercepts=qlogis(seq(.9,.1,by=-.1)),
                    values=1:10, interceptRef=1, cumprob=plogis, inverse=qlogis,
                    conf.int=0, method='interpolated')
for(a in c(.01, seq(0, 1, by=.05), .99))
  cat(a, qu(a, qlogis(.9)), '\n')


set.seed(3)
n <- 300
x1 <- runif(n)
ddist <- datadist(x1); options(datadist='ddist')
y  <- x1 + runif(n)
x1[1:35] <- NA
dat <- data.frame(x1, y)
f <- orm(y ~ x1, x=TRUE, y=TRUE)
f
g <- bootcov(f, B=50)
x1s <- seq(0, 1, by=.1)
pg  <- Predict(g, x1=x1s, boot.type='basic')
cof <- c(coef(f)[f$interceptRef], coef(f)[length(coef(f))])
cof
apply(g$boot.Coef, 2, mean)
sqrt(var(g$boot.Coef[,2]))
a <- aregImpute(~ x1 + y, data=dat)
h <- fit.mult.impute(y ~ x1, orm, a, data=dat)
pmi <- Predict(h, x1=x1s)

plot(Predict(f, x1=x1s),
             addpanel=function(...) {
               with(pg, {llines(x1,  lower, col='red')
                         llines(x1,  upper, col='red')
                         lpoints(x1, yhat,  col='red')
                         llines(x1s,  cof[1] + cof[2]*x1s, col='green')
                       })
               with(pmi, {lpoints(x1, lower, col='black')
                          lpoints(x1, upper, col='black')
                          lpoints(x1, yhat,  col='black')})
             })


require(rms)
getHdata(nhgh)
w <- subset(nhgh, age >= 21 & dx==0 & tx==0,
            select=-c(dx,tx))
dd <- datadist(w); options(datadist='dd')
set.seed(1)
w$ghe <- w$gh + round(runif(length(w$gh), -.5, .5), 2)

Ecdf(~ gh, groups=is.na(sub), data=w)
Ecdf(~ age, groups=is.na(sub), data=w)
with(w, table(is.na(sub), gh))
length(unique(w$gh))
with(w, tapply(gh, is.na(sub), function(x) length(unique(x))))
w2 <- subset(w, !is.na(sub))
wdata <- 2
## If substitute ghr for gh, boot problem goes away for w2
g <- orm(gh ~ age, family='loglog', data=if(wdata == 1) w else w2,
         x=TRUE, y=TRUE)
set.seed(2)
gb <- bootcov(g, B=100, pr=TRUE)
ages <- seq(25, 80, by=5)
bootcl  <- Predict(gb, age=ages, boot.type=c('percentile','basic')[2])
bootclcov <- Predict(gb, age=ages, usebootcoef=FALSE)
X <- predict(gb, newdata=bootcl, type='x')
br <- gb$boot.Coef[,1] + X %*% t(gb$boot.Coef[,-1])
if(wdata == 1) br1 <- br else br2 <- br
z <- quantile(br[1,], c(.025,.975))
plot(Predict(g, age=ages), ylim=c(-1.5,1.5), addpanel=function(...) {
  lpoints(23, z, col='red')
  for(j in 1:12) lpoints(ages[j], br[j,], col=gray(.9))
  with(bootcl, {llines(age, lower, col='blue')
                llines(age, upper, col='blue')
                lpoints(age, yhat, col='blue')})
  with(bootclcov, {llines(age, lower, col='red')
                   llines(age, upper, col='red')})
})

# For age of 70 "manually" find predicted median
# was predict(f, ...)  why?
p <- predict(g, newdata=data.frame(age=70), type='fitted.ind')
cumsum(p)
median(rep(g$yunique, round(1000000*p)))  # 5.6
xb <- Function(g)(age=70)
intercepts <- coef(f)[1 : num.intercepts(f)]
# Compute Prob(Y <= y) from Prob(Y >= y) by shifting one level
# Prob(Y > y) = Prob(Y >= y + epsilon)
cumprob <- f$trans$cumprob
# xless(cumprob(intercepts + xb))
names(intercepts) <- Lag(names(intercepts))
names(intercepts) <- gsub('>=', '>', names(intercepts))
intercepts
probYley <- 1 - cumprob(intercepts + xb)
names(probYley) <- gsub('>', '<=', names(probYley))
probYley   # 5.6 gives prob Y <= 5.6 = .50899.  Interpolated median 5.59

# pgty <- f$trans$cumprob(intercepts + xb)
# Prob(Y <= y) = Prob(Y < y + epsilon) = 1 - Prob(Y >= y + epsilon)
pleq <- cumprob(coef(f)[1:num.intercepts(f)] + xb)
lp <- coef(f)[f$interceptRef] + xb

## Look at bootstrap variation in median gh for both subsets
B <- 2000; meds1 <- meds2 <- numeric(B)
y1 <- w$gh;    y2 <- w2$gh
n1 <- nrow(w); n2 <- nrow(w2)
pb <- setPb(B, every=50)
for(i in 1:B) {
  pb(i)
  s <- sample(1:n1, n1, replace=TRUE)
  meds1[i] <- median(y1[s])
  s <- sample(1:n2, n2, replace=TRUE)
  meds2[i] <- median(y2[s])
}
table(meds1); table(meds2)


# See how to check intercepts against linear model assumptions
require(rms)
set.seed(1)
n <- 1000
x1 <- runif(n)
y <- 30 + x1 + rnorm(n)
f <- orm(y ~ x1, family=probit)
y2 <- y + 20
f2 <- orm(y2 ~ x1, family=probit)
plot(coef(f), coef(f2)) # unaffected by shift

g  <- ols(y ~ x1)
yu <- f$yunique[-1]
ns <- num.intercepts(f)
s  <- g$stats['Sigma']
alphas <- coef(f)[1:ns]
plot(-yu/s, alphas, type='l',
     xlab=expression(-y/s), ylab=expression(alpha[y]))
co <- coef(lm.fit(cbind(1, -yu/s), alphas))
text(-32, 2, paste('Slope:', round(co[2], 4)))
abline(a=co[1], b=co[2], col='gray70')

## Compare coefficients with those from partial likelihood (Cox model)
orm(y ~ pol(x1,2), family=loglog)
cph(Surv(y) ~ pol(x1,2))


## Simulate from a linear model with normal residuals and compute
## quantiles for one x value, two ways

set.seed(7)
n <- 10000
x <- rnorm(n)
y <- round(x + rnorm(n), 2)
f <- ols(y ~ x)
k <- coef(f)
s <- f$stats['Sigma']
print(c(k, s))
k[1] + qnorm((1:3)/4) * s

g <- orm(y ~ x, family='probit')
quant <- Quantile(g)
lp <- predict(g, data.frame(x=0))
for(qu in (1:3)/4) print(quant(qu, lp))
