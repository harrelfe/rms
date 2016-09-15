require(rms)
#
# Time to progression/death for patients with monoclonal gammopathy
# Competing risk curves (cumulative incidence)
# status variable must be a factor with first level denoting right censoring
m <- upData(mgus1, stop = stop / 365.25, units=c(stop='years'),
            labels=c(stop='Follow-up Time'), subset=start == 0)
f <- npsurv(Surv(stop, event) ~ sex, data=m)
levels(m$event)
f$n.risk
## Compute long-hand for checking
times <- f$time
sex <- c(rep('female', f$strata['sex=female']),
         rep('male',   f$strata['sex=male']))
n <- length(times)
nrisk <- numeric(n)
fu <- m$stop
for(i in 1 : n) nrisk[i] <- sum(fu >= times[i] - 1e-7 & m$sex == sex[i])
w <- data.frame(sex, times, nrisk, f$n.risk)
w
## xless(w)


times <- seq(0, 36, by=4)
g <- summary(f, times=times, print.it=FALSE)
## unclass(g)
sex <- as.character(g$strata)
n <- length(times)
nrisk <- matrix(0, nrow=length(times), ncol=2)
colnames(nrisk) <- c('female', 'male')
for(sx in c('female', 'male')) 
  for(i in 1 : n) nrisk[i, sx] <- sum(fu >= times[i] - 1e-7 & m$sex == sx)
nrisk


par(mar=c(8, 4, 1, 1))
survplot(f, state='pcm', n.risk=TRUE, xlim=c(0, 20), ylim=c(0, .5), col=1:2,
         y.n.risk=-.15)

survplotp(f, state='pcm', xlim=c(0, 20), ylim=c(0, .5))



survplot(f, state='death', add=TRUE, col=3)

f <- npsurv(Surv(stop, event) ~ sex, data=m)
survplot(f, state='death', n.risk=TRUE, conf='diffbands')
