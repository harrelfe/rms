## See http://stats.stackexchange.com/questions/147733/equation-of-a-fitted-smooth-spline-and-its-analytical-derivative

## This fits a natural spline (linear tail restricted) using the truncated
## power basis.  Default knots are not used; instead specify 4 knots.
require(rms)
x <- 1:11
y <- c(0.2,0.40, 0.6, 0.75, 0.88, 0.99, 1.1, 1.15, 1.16, 1.16, 1.16 )
dd <- datadist(x); options(datadist='dd')

f <- ols(y ~ rcs(x, c(3, 5, 7, 9)))
f
anova(f)
ggplot(Predict(f)) + geom_point(aes(x=x, y=y), data=data.frame(x,y))
Function(f)   ## if have latex installed can also use latex(f)

## Function re-expresses the restricted cubic spline in simplest form

## The first derivative is:

## function(x) 0.174 - 3 * 0.00279 * pmax(x - 3, 0) ^ 2 + 3 * 0.0015 * pmax(x - 5, 0) ^ 2 + ...
