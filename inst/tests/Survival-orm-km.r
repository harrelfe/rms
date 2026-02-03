# For many simulated discrete survival and censoring times make sure
# that survival::survfit and Survival.orm agree with respect to distinct failure times
# and survival probabilities at the observed failure times and at
# a regular grid of times including times halfway between observed times.
# Also check agreement of confidence limits when survfit is told to use log-log basis.

require(rms)
require(survival)

# Simple example with no early censoring, with last time occurring as
# both censored and uncensored
y <- c(1:4, 4)
e <- c(1, 1, 1, 1, 0)
Y <- Ocens(y, ifelse(e == 1, y, Inf))
w <- Ocens2ord(Y)
w
unclass(w)
S <- Surv(y, e)
summary(survfit(S ~ 1, conf.type='log-log'))
f <- orm.fit(y=Y, trace=1, family='loglog')
f$yunique
f$coefficients
exp(-exp(-f$coefficients))
s <- Survival(f)
s()

nignore <- 0
set.seed(21)
pr  <- function(...) {}
chk <- function(a, b, msg) if(! isTRUE(all.equal(a, b))) stop(msg)

for(i in 1: 2500) {
  pr('i=', i, '')
  n    <- sample(10:40, 1)
  t    <- sample(1 : n, n, replace=TRUE)
  cens <- sample(1 : n, n, replace=TRUE)
  y    <- pmin(t, cens)
  e    <- ifelse(t <= cens, 1, 0)
  # If there are any right-censored points before all uncensored times, these
  # are set to NA because they don't contribution to the likelihood
  # Skip these cases
  if(any(y[e == 0] < min(y[e == 1]))) {
    nignore <- nignore + 1
    next
  }
  Y    <- Ocens(y, ifelse(e == 1, y, Inf))
  S    <- Surv(y, e)
  f    <- orm.fit(y=Y, family='loglog')
  s    <- Survival(f)
  pr('a')
  s1   <- s(conf.int=0.95)
  t1   <- s1$time
  low1 <- s1$lower
  up1  <- s1$upper
  s1   <- s1$surv
  g   <- survfit(S ~ 1, conf.type='log-log')
  u   <- summary(g)
  t2  <- u$time
  s2  <- u$surv
  low2 <- u$lower
  up2  <- u$upper
  stopifnot(all(is.na(low2) == is.na(up2)))

  chk(t1, t2, 'unequal times')
  chk(s1, s2, 'unequal survival probs')
  m <- max(abs(f$u))
  stopifnot(m < 1e-13)

  chk(low1, low2, 'lower CLs not equal')
  chk(up1,  up2,  'upper CLs not equal')

  ts <- seq(0, max(y) + 2, by=0.5)
  pr('b\n')
  s1 <- s(times=ts, conf.int=0.95)
  low1 <- s1$lower
  up1  <- s1$upper
  s1   <- s1$surv
  s2 <- summary(g, times=ts, extend=TRUE)

  low2 <- s2$lower
  up2  <- s2$upper
  s2   <- s2$surv

  # survfit erroneously extends non-zero S(t) beyond observed points if last censored
  # point is at or beyond the highest uncensored point; set to NA like Survival.orm
  # Note: If the highest point has an uncensored point at that value, S(t >= highest) = 0
  maxc <- max(y[e == 0])
  maxu <- max(y[e == 1])
  if(any(e == 0) && maxc >= maxu) {
    s2  [ts > maxc] <- NA
    low2[ts > maxc] <- NA
    up2 [ts > maxc] <- NA
  }
  chk(s1,     s2, 'unequal survival probs at requested times')
  chk(low1, low2, 'unequal lower CLs at requested times')
  chk(up1,  up2,  'unequal upper CLs at requested times')
}

cat(nignore, 'iterations ignored\n')


