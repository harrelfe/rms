require(rms)
n <- 1000
cvd <- data.frame(id    = sample(1 : 100, n, TRUE),
                  bmi   = rnorm(n, 25, 3),
                  wtpre = sample(c(-1, 1), n, TRUE))
t.enter <- runif(n)
t.exit  <- t.enter + runif(n)
cens    <- sample(0:1, n, TRUE)
cvd$S <- Surv(t.enter, t.exit, cens)

label(cvd$id) <- "Id"
label(cvd$bmi) <- "BMI"
units(cvd$bmi) <- "Kg/m2"

dd <- datadist(cvd); options(datadist = "dd")
#S <- with(cvd, Surv(t.enter, t.exit, cens))

cph(S ~ rcs(bmi) + cluster(id) + wtpre, data=cvd)
f <- cph  (S ~ pol(bmi, 3) + cluster(id) + offset(wtpre), data = cvd, eps=1e-6)
g <- coxph(S ~ pol(bmi, 3) + cluster(id) + offset(wtpre), data = cvd,
           control=coxph.control(eps=1e-6))
coef(f) - coef(g)
f
g

d <- data.frame(bmi=20, wtpre=5)
Predict(f, bmi=20)  # will not work
predict(f, d)       # 5.031003
k <- coef(f)
with(d, - f$center + k[1]*bmi + k[2] * bmi^2 + k[3] * bmi^3 + wtpre)  # 5.031003
predict(g, d)   # bombs

f <- cph(S ~ pol(bmi, 3) + offset(wtpre), data=cvd, eps=1e-6)
g <- coxph(S ~ pol(bmi, 3) + offset(wtpre), data=cvd)
coef(f) - coef(g)
predict(f, d)

z <- expand.grid(x1=LETTERS[1:3], x2=c('a','b','c'), reps=1:10)
z$S <- Surv(runif(nrow(z)))
cph(S ~ x1 * strat(x2), data=z)
