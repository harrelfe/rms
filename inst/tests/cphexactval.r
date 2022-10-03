require(rms)
data(cancer, package='survival')
f <- cph(Surv(time, status) ~ age, data=cancer, x=TRUE, y=TRUE, surv=TRUE,
         method='exact')
validate(f, B=20)
