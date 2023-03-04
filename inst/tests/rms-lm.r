# From Shawn Garbett 2022-11-21
require(rms)
x <- (1:200 + rnorm(200)) / 50
y <- 5*sin(x) + x^2 + rnorm(100)
plot(x, y, main="test data")
lines(x, 5*sin(x)+x^2)
model <- lm(y ~ rcs(x,3))
summary(model)
if(require(strip)) {
  stripped <- strip(model, "predict")
  x <- seq(0, 4, by=0.2)
  plot(x, predict(stripped, data.frame(x=x)), main="Predicted vs. Truth",
       ylab="Predicted Y")
  lines(x, 5*sin(x)+x^2, col='red')
}

