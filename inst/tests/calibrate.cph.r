## Eduwin Pakpahan <epakpahan@gmail.com>

require(rms)
require(haven)

d <- read_dta("pakpahan.dta")

fit <- cph(Surv(data_dftime, data_demfu) ~ data_age, ties="breslow", data=d,
           surv=TRUE, x=T, y=T, time.inc=1200)
fit

cal <- calibrate(fit, u=1200, B=120)
plot(cal, subtitles=FALSE)

cal_KM <- calibrate(fit, cmethod='KM', u=1200, m=10, B=40)

plot(cal_KM, add=TRUE)
