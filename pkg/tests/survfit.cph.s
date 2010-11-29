plots      <- FALSE
testDesign <- FALSE

## Simulate a small example to compare results with survival package
n <- 100
set.seed(1)
age <- 50 + 12*rnorm(n)
sex <- factor(sample(c('Male','Female'), n, rep=TRUE))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
e <- ifelse(dt <= cens, 1, 0)
dt <- round(pmin(dt, cens))
dtmax <- tapply(dt, sex, max)
require(survival)
Srv <- Surv(dt, e)
f <- coxph(Srv ~ age*strata(sex))
levels(sex)
new <- expand.grid(age=50, sex=levels(sex))

new1 <- new[1,]
new2 <- new[2,]

if(plots)
  {
    gr <- function(col=gray(.9), lwd=1)
      abline(h=seq(.2,.8,by=.2), col=col, lwd=lwd)
    par(mfrow=c(2,2))

    s <- survfit(f, new1)
    plot(s, conf.int=TRUE, main=new1$sex)
    lines(s, col='red')
    gr()
    s <- survfit(f, new2)
    plot(s, conf.int=TRUE, main=new2$sex)
    lines(s, col='red')
    gr()
  }

s  <- survfit(f, new)
st <- rep(names(s$strata), s$strata)
i <- 0
for(sx in levels(sex))
  {
    i <- i + 1
    cat(sx, '\t', 'survfit.coxph\n')
    j <- st==paste('sex=', sx, sep='')
    z <- with(s, data.frame(time=time[j], surv=surv[j,i],
                            std.err=std.err[j,i],
                            lower=lower[j,i],
                            upper=upper[j,i]))
    print(z)
  }

options(digits=3)
z <- with(s, data.frame(strata=st, time=time, surv=surv, se=std.err,
                        lower=lower, upper=upper))
if(plots)
  {
    with(subset(z, strata=='sex=Female'),
         {
           plot(c(time,dtmax['Female']), c(surv.1, min(surv.1)),
                type='s', col='red', xlim=c(0,15), ylim=c(0,1))
           points(time, surv.1, pch='f', col='red')
           lines(time, lower.1, type='s', col='red')
           lines(time, upper.1, type='s', col='red')
         })
    with(subset(z, strata=='sex=Male'),
         {
           lines(c(time,dtmax['Male']), c(surv.2, min(surv.2)),
                 type='s', col='green')
           points(time, surv.2, pch='m', col='green')
           lines(time, lower.2, type='s', col='green')
           lines(time, upper.2, type='s', col='green')
         })
    gr()
  }

require(rms)
dd <- datadist(age,sex); options(datadist='dd')
Srv <- Surv(dt, e)
g <- cph(Srv ~ age*strat(sex), surv=TRUE)

for(sx in levels(sex))
  {
    k <- survest(g, data.frame(age=50, sex=sx))
    cat(sx, '\t', 'survest surv=TRUE\n')
    print(with(k, data.frame(time=time, surv=surv,
                             std.err=std.err, lower=lower, upper=upper)))
  }

if(plots)
  {
    survplot(g, sex=., age=50, conf.int=TRUE)
    w <-survest(g, data.frame(age=50, sex='Female'))
    points(w$time, w$surv, pch='f')
    w <- survest(g, data.frame(age=50, sex='Male'))
    points(w$time, w$surv, pch='m')
    gr(col='blue', lwd=.2)
  }

h <- cph(Srv ~ age*strat(sex), x=TRUE, y=TRUE)
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
    survplot(h, sex=., age=50, conf.int=TRUE)
    w <- survest(h, data.frame(age=50, sex='Female'))
    points(w$time, w$surv, pch='f')
    w <- survest(h, data.frame(age=50, sex='Male'))
    points(w$time, w$surv, pch='m')
    gr(col='blue', lwd=.2)
  }

if(testDesign)
  {
    ## To compare with Design
    require(Design)
    Srv <- Surv(dt, e)
    new <- expand.grid(age=50, sex=levels(sex))
    dd <- datadist(age,sex); options(datadist='dd')
    ## Temporary fix until Design updated around 10Sep09
    Varcov <- function(object, ...) UseMethod("Varcov")
    Varcov.cph <- function(object, regcoef.only=FALSE, ...)
      Varcov.default(object, regcoef.only, ...)
    Varcov.default <- function(fit, regcoef.only=FALSE)
      {
        vc <- fit$Varcov
        if(length(vc)) {
          if(regcoef.only) return(fit$var) else
          return(vc(fit,which='var'))
        }
        cov <- fit$var
        if(is.null(cov)) stop("fit does not have variance-covariance matrix")
        if(regcoef.only)
          {
            p <- length(fit$coefficients)  # 14Sep00
            cov <- cov[1:p, 1:p, drop=FALSE]
          }
        cov
      }
  
    g <- cph(Srv ~ age*strat(sex), surv=TRUE)

    options(digits=3)
    for(sx in levels(sex))
      {
        k <- survest(g, data.frame(age=50, sex=sx))
        cat(sx, '\t', 'survest Design surv=TRUE\n')
        print(with(k, data.frame(time=time, surv=surv,
                                 std.err=std.err, lower=lower, upper=upper)))
      }
    g <- cph(Srv ~ age*strat(sex), x=TRUE, y=TRUE)
    cat('cph x=T y=T survfit\n')
    print(unclass(survfit(g, new, conf.type='log')))
    
  }
