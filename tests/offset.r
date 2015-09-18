library(rms)

cvd <- read.csv("~/tmp/cvd.csv")

label(cvd$id) <- "Id"
label(cvd$bmi) <- "BMI"
units(cvd$bmi) <- "Kg/m2"

dd <- datadist(cvd); options(datadist = "dd")
S <- with(cvd, Surv(t.enter, t.exit, cens))

tt <- function() {source('~/R/rms/R/rms.s');source('~/R/rms/R/cph.s');source('~/R/rms/R/predictrms.s');source('~/R/rms/R/rmsMisc.s');source('~/R/rms/R/rms.trans.s')}

cph(S ~ rcs(bmi) + cluster(id) + wtpre, data=cvd)
f <- cph  (S ~ pol(bmi, 3) + cluster(id) + offset(wtpre), data = cvd, eps=1e-6)
g <- coxph(S ~ pol(bmi, 3) + cluster(id) + offset(wtpre), data = cvd,
           control=coxph.control(eps=1e-6))
coef(f) - coef(g)
f
g

d <- data.frame(bmi=20, wtpre=5) #, id=3000619)
Predict(f, bmi=20)  # will not work
predict(f, d)       # 5.1806
k <- coef(f)
with(d, - f$center + k[1]*bmi + k[2] * bmi^2 + k[3] * bmi^3 + wtpre)  # 5.1806
predict(g, d)   # bombs

f <- cph(S ~ pol(bmi, 3) + offset(wtpre), data=cvd, eps=1e-6)
g <- coxph(S ~ pol(bmi, 3) + offset(wtpre), data=cvd)
coef(f) - coef(g)
predict(f, d)

z <- expand.grid(x1=LETTERS[1:3], x2=c('a','b','c'), reps=1:10)
z$S <- Surv(runif(nrow(z)))
cph(S ~ x1 * strat(x2), data=z)
