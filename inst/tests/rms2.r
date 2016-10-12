## From Jerome Asselin  https://github.com/harrelfe/rms/issues/32
require(rms)
df <- data.frame(y=rnorm(21), day=weekdays(Sys.Date()-(1:21), abbr=TRUE))
df$day.ordered <-
  with(df, factor(as.character(day),
                  levels=c("Sun","Mon","Tue","Wed","Thu","Fri","Sat"),
                  ordered=TRUE))
options(contrasts=c("contr.treatment", "contr.treatment"))
fit1 <- ols(y ~ day, data=df)
fit2 <- ols(y ~ day.ordered, data=df)
df.char <- df
df.char$day.ordered <- as.character(df.char$day.ordered)
w <- cbind(orig          = predict(fit1),
      orig.newdata       = predict(fit1, newdata=df),
      ordered            = predict(fit2),
      ordered.newdata    = predict(fit2, newdata=df),
      ordered.workaround = predict(fit2, newdata=df.char))
round(w[, -1] - w[, 1], 3)
