require(rms)
require(survival)
data(cancer)
f <- cph(Surv(time, status) ~ age, data=cancer, x=TRUE, y=TRUE, surv=TRUE,
         method='exact')
validate(f, B=20)
