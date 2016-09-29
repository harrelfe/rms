# Survival time with stratification.  Thanks: Cathy Jenkins
require(rms)
Load(sampledf)
S  <- with(sampledf, Surv(fu, death))
dd <- datadist(sampledf); options(datadist='dd')

f <- cph(S ~ rcs(age, c(1, 2, 4)) + rcs(sqrt(cd4), sqrt(c(210,475,875))) +
           strat(site), data=sampledf, x=TRUE, y=TRUE, iter.max=30, eps=1e-10)
g <- coxph(S ~ rcs(age, c(1, 2, 4)) + rcs(sqrt(cd4), sqrt(c(210,475,875))) +
             strata(site), data=sampledf, x=TRUE, y=TRUE, control=coxph.control(eps=1e-10, iter.max=30))

# -------------------------- #
#   Survival probabilities   #
#         at 1 year          #
#  for fixed ages/CD4 counts #
# -------------------------- #
Pd <- expand.grid(age=c(1,2,6,9,15),
                  cd4=c(100,200,350,500),
                  site=levels(sampledf$site))
pd <- Pd[1, ]

a <- survfit(f, newdata=pd); a$strata
b <- survfit(g, newdata=pd); b$strata

h <- function(a, b, chkeq=FALSE) {
  a <- sapply(a, length)
  b <- sapply(b, length)
  k <- unique(c(names(a), names(b)))
  z <- matrix(NA, nrow=length(k), ncol=2, dimnames=list(k, c('a', 'b')))
  z[names(a), 'a'] <- a
  z[names(b), 'b'] <- b
  print(z)
  if(chkeq) {
    k <- intersect(names(a), names(b))
    for(n in k)
      cat(n, ' equal:', all.equal(a[[n]], b[[n]]), '\n')
    }
}
h(a, b)

a <- survfit(f, newdata=Pd)
b <- survfit(g, newdata=Pd)
h(a, b, chkeq=TRUE)

z <- summary(survfit(g, newdata=Pd[1:1,]), times=5)


z <- survest(f, newdata=Pd[33,], times=5)

comp <- function(a, b, ntimes=1, time=1, ib=TRUE) {
  b$std.err <- b$std.err / b$surv
  for(n in c('time', 'surv', 'std.err', 'lower', 'upper',
             if(length(a$strata)) 'strata')) {
    x <- a[[n]]
    y <- b[[n]][ib]
    if(n %nin% c('time', 'strata') && ntimes > 1) {
      x <- x[, time]
      y <- y[seq(time, length(y), by=ntimes)]
      }
    cat(n, ' equal:',
        if(length(x) == length(y)) all.equal(x, y)
        else
          paste('lengths:', length(x), length(y)),
        '\n')
  }
}

chk <- function(f, g, strat=FALSE) {
  a <- survest(f, newdata=Pd[33,], times=5)
  b <- summary(survfit(g, newdata=Pd[33,]), times=5)
  cat('-------------------------- newdata 1 row, 1 time\n')
  comp(a, b, ib=if(strat) 2 else TRUE)

  a <- survest(f, newdata=Pd, times=5)
  b <- summary(survfit(g, newdata=Pd), times=5)
  cat('-------------------------- newdata all, 1 time\n')
  comp(a, b)

  a <- survest(f, newdata=Pd, times=5:6)
  b <- summary(survfit(g, newdata=Pd), times=5:6)
  cat('-------------------------- newdata all, 2 times\n')
  comp(a, b, ntimes=2, time=1)
}
chk(f, g, strat=TRUE)

## Try with no strata

f <- cph  (S ~ rcs(age, c(1, 2, 4)) + rcs(sqrt(cd4), sqrt(c(210,475,875))) +
           site, data=sampledf, x=TRUE, y=TRUE, iter.max=30, eps=1e-10)
g <- coxph(S ~ rcs(age, c(1, 2, 4)) + rcs(sqrt(cd4), sqrt(c(210,475,875))) +
             site, data=sampledf, x=TRUE, y=TRUE, control=coxph.control(eps=1e-10, iter.max=30))
cbind(coef(f), coef(g))

chk(f, g)
