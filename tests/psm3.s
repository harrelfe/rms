# From IM Nolte <i.m.nolte@umcg.nl>
require(rms)
set.seed(1)
n <- 1000
v    <- rbinom(n, 2, 0.2)
time <- rnorm(n, v / 10 + 2, 0.5)
c <- ifelse(time < 0.5, 2, ifelse(time > 3.5, 0, ifelse(time > 2.5, 3, 1)))
time[c==2] <- 0.5
time[c==0] <- 3.5
time2 <- time + 0.1
time2[c==3] <- time[c==3] + runif(sum(c==3), 0.1, 0.5)
S <- Surv(time, time2, c, type="interval")
survreg(S ~ v, dist='gaussian')
psm(S ~ v, dist='gaussian')
