require(rms)
set.seed(1)
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

Predict(f, bmi=20, offset=list(wtpre=5))   # 4.849534
d <- data.frame(bmi=20, wtpre=5)
predict(f, d)       # 4.849534
k <- coef(f)
mp <- function(fit, bmi=20, wtpre=5) {
  k  <- coef(fit)
  k0 <- if(length(fit$center)) - fit$center
  else {
    m  <- fit$means
    k0 <- - unname(k[1] * m[1] + k[2] * m[2] + k[3] * m[3])
    }
  k0 + k[1]*bmi + k[2] * bmi^2 + k[3] * bmi^3 + wtpre
}
mp(f)   # 4.849534
mp(g)   # "

f <- cph(S ~ pol(bmi, 3) + offset(wtpre), data=cvd, eps=1e-6)
g <- coxph(S ~ pol(bmi, 3) + offset(wtpre), data=cvd,
           control=coxph.control(eps=1e-6))
coef(f) - coef(g)
mp(f)   # 4.849534
mp(g)   # "
predict(f, d)   # -.1504; ignores offset
Predict(f, bmi=20, offset=list(wtpre=5))   # -.1504

p1 <- Predict(f, bmi, offset=list(wtpre=0))
p2 <- Predict(f, bmi, offset=list(wtpre=5))
plot(rbind(p1, p2))

z <- expand.grid(x1=LETTERS[1:3], x2=c('a','b','c'), reps=1:10)
z$S <- Surv(runif(nrow(z)))
cph(S ~ x1 * strat(x2), data=z)
