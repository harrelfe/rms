## Eduwin Pakpahan <epakpahan@gmail.com>

if(require(haven)) {
  require(rms)

  d <- read_dta("pakpahan.dta")

  fit <- cph(Surv(data_dftime, data_demfu) ~ data_age, method="breslow", data=d,
             surv=TRUE, x=T, y=T, time.inc=1200)
  print(fit)

  cal <- calibrate(fit, u=1200, B=120)
  plot(cal, subtitles=FALSE)

  cal_KM <- calibrate(fit, cmethod='KM', u=1200, m=10, B=40)

  plot(cal_KM, add=TRUE)
}
