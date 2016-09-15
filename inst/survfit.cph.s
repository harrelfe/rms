## Compare SE of log survival probability from survest and survfit

require(rms)
set.seed(123)
n <- 200
age <- rnorm(n, 50, 10)
x <- 50*runif(n)
ct <- round(365*runif(n))
h <- .003*exp(.005*age+0.008*x)
ft <- round(-log(runif(n))/h)
status <- ifelse(ft <= ct,1,0)
ft <- pmin(ft, ct)

S <- Surv(ft, status)
fit <- cph(S ~ age + x, x=TRUE, y=TRUE)
d <- data.frame(age=mean(age), x=mean(x))
s <- survest(fit, d, times=56)
prn(with(s, cbind(time, surv, std.err, lower, upper)), 'survest')
s <- survfit(fit, d)

k <- which(s$time == 56)
prn(with(s, cbind(time, surv, std.err, lower, upper)[k,]), 'survfit')

fit <- cph(S ~ age + x, x=TRUE, y=TRUE, surv=TRUE, time.inc=56)
k <- which(fit$time == 56)
prn(fit$surv.summary, 'cph surv=T')
s <- survest(fit, d, times=56)
prn(with(s, cbind(time, surv, std.err, lower, upper)),
    'survest from cph surv=T')
s <- survfit(fit, d)
k <- which(s$time == 56)
prn(with(s, cbind(time, surv, std.err, lower, upper)[k,]),
    'survfit from cph surv=T')


survest(fit, data.frame(age=40, x=25), times=56)
 
pp <- rms:::survfit.cph(fit, data.frame(age=40, x=25),se.fit=TRUE)
cbind(pp$std.err, pp$lower,pp$upper)[pp$time==56]
 

##--------------------------------------------------------------

require(survival)

plots      <- TRUE
topdf      <- TRUE
testrms    <- FALSE
testDesign <- FALSE
additive   <- FALSE
roundt     <- FALSE
chkpts     <- FALSE

## Simulate a small example to compare results with survival package
nfemale <- 100
nmale   <- 9*nfemale
n <- nfemale + nmale
set.seed(1)
age <- 50 + 12*rnorm(n)
sex <- as.factor(c(rep('Male', nmale), rep('Female', nfemale)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+1.5*(sex=='Female'))
dt <- -log(runif(n))/h
e <- ifelse(dt <= cens, 1, 0)
dt <- if(roundt) round(pmin(dt, cens)) else pmin(dt, cens)
dtmax <- tapply(dt, sex, max)
Srv <- Surv(dt, e)
f <- coxph(if(additive) Srv ~ age + strata(sex) else Srv ~ age*strata(sex))
levels(sex)
new <- expand.grid(age=50, sex=levels(sex))
new1 <- new[1,]
new2 <- new[2,]

if(plots)
  {
    if(topdf) pdf('/tmp/z.pdf')
    gr <- function(col=gray(.9), lwd=1)
      abline(h=seq(.2,.8,by=.2), col=col, lwd=lwd)
    par(mfrow=c(2,2))

    s <- survfit(f, new1, censor=FALSE)
    plot(s, conf.int=TRUE,
         main=paste(new1$sex, 'coxph survfit newdata new1'))
    lines(s, col='red')
    gr()
    s <- survfit(f, new2, censor=FALSE)
    plot(s, conf.int=TRUE,
         main=paste(new2$sex, 'coxph survfit newdata new2'))
    lines(s, col='red')
    gr()
  }

s  <- survfit(f, new, censor=FALSE)
plot(s, main='coxph combined newdata plot.survfit')
gr()

z <- with(s, data.frame(time, surv, std.err, lower, upper, se=std.err,
                        strata=c(rep('Female', s$strata[1]),
                                 rep('Male',   s$strata[2]))) )
z

options(digits=3)

if(FALSE && plots)
  {
    with(subset(z, strata=='Female'),
         {
           plot(c(time,dtmax['Female']), c(surv.1, min(surv.1)),
                type='s', col='red', xlim=c(0,15), ylim=c(0,1),
                xlab='', ylab='')
           if(chkpts) points(time, surv.1, pch='f', col='red')
           lines(time, lower.1, type='s', col='red')
           lines(time, upper.1, type='s', col='red')
         })
    with(subset(z, strata=='Male'),
         {
           lines(c(time,dtmax['Male']), c(surv.2, min(surv.2)),
                 type='s', col='green')
           if(chkpts) points(time, surv.2, pch='m', col='green')
           lines(time, lower.2, type='s', col='green')
           lines(time, upper.2, type='s', col='green')
         })
    title('coxph combined newdata manual')
    gr()
  }

if(testrms)
  {
    require(rms)
    system('cat ~/R/rms/pkg/R/*.s > /tmp/rms.s')
    source('/tmp/rms.s')
    dd <- datadist(age,sex); options(datadist='dd')
    Srv <- Surv(dt, e)
    g <- cph(if(additive) Srv ~ age + strat(sex) else Srv ~ age*strat(sex),
             surv=TRUE)
    
    for(sx in levels(sex))
      {
        k <- survest(g, data.frame(age=50, sex=sx))
        cat(sx, '\t', 'survest surv=TRUE\n')
        print(with(k, data.frame(time=time, surv=surv,
                                 std.err=std.err, lower=lower, upper=upper)))
      }
    
    if(plots)
      {
        survplot(g, sex, age=50, conf.int=TRUE)
        w <-survest(g, data.frame(age=50, sex='Female'))
        if(chkpts) points(w$time, w$surv, pch='f')
        w <- survest(g, data.frame(age=50, sex='Male'))
        if(chkpts) points(w$time, w$surv, pch='m')
        title('rms survplot + survest surv=T')
        gr()
      }
    
    h <- cph(if(additive) Srv ~ age + strat(sex) else Srv ~ age*strat(sex),
             x=TRUE, y=TRUE)
    s <- survfit(h, new)
    unclass(s)
    st <- rep(names(s$strata), s$strata)
    i <- 0
    for(sx in levels(sex))
      {
        i <- i + 1
        cat(sx, '\t', 'survfit.cph surv=F\n')
        j <- st==paste('sex=', sx, sep='')
        z <- with(s, data.frame(time=time[j], surv=surv[j,i],
                                std.err=std.err[j,i],
                                lower=lower[j,i],
                                upper=upper[j,i]))
        print(z)
        k <- survest(h, data.frame(age=50, sex=sx))
        cat(sx, '\t', 'survest surv=F\n')
        print(with(k, data.frame(time=time, surv=surv,
                                 std.err=std.err, lower=lower, upper=upper)))
      }

## i <- s2$strata
## with(s2, data.frame(time=time[,i], surv=surv[,i], se=std.err[,i],
##                          lower=lower[,i], upper=upper[,i]))

    if(plots)
      {
        survplot(h, sex, age=50, conf.int=TRUE)
        w <- survest(h, data.frame(age=50, sex='Female'))
        if(chkpts) points(w$time, w$surv, pch='f')
        w <- survest(h, data.frame(age=50, sex='Male'))
        if(chkpts) points(w$time, w$surv, pch='m')
        title('rms survplot + survest x=T y=T')
        gr()
      }
  }

if(testDesign)
  {
    ## To compare with Design
    require(Design)
    Srv <- Surv(dt, e)
    new <- expand.grid(age=50, sex=levels(sex))
    dd <- datadist(age,sex); options(datadist='dd')
  
    g <- cph(if(additive) Srv ~ age + strat(sex) else Srv ~ age*strat(sex),
             surv=TRUE)
    if(plots)
      {
        survplot(g, sex=NA, age=50, conf.int=TRUE, conf='bands')
        gr()
        title('Design survplot surv=T')
      }

    options(digits=3)
    for(sx in levels(sex))
      {
        k <- survest(g, data.frame(age=50, sex=sx))
        cat(sx, '\t', 'survest Design surv=TRUE\n')
        print(with(k, data.frame(time=time, surv=surv,
                                 std.err=std.err, lower=lower, upper=upper)))
      }
    g <- cph(if(additive) Srv ~ age + strat(sex) else Srv ~ age*strat(sex),
             x=TRUE, y=TRUE)
    cat('cph x=T y=T survfit\n')
    print(unclass(survfit(g, new, conf.type='log')))
    if(plots)
      {
        survplot(g, sex=NA, age=50, conf.int=TRUE, conf='bands')
        title('Design survplot x=T y=T')
        gr()
      }
  }
if(topdf) dev.off()
