# From http://stats.stackexchange.com/questions/195198
require(rms)
d1 <- data.frame(cohort='one', sex='male', y=c(.476,
.84,
1.419,
0.4295,
0.083,
2.9595,
4.20125,
1.6605,
3.493,
5.57225,
0.076,
3.4585))
d2 <- data.frame(cohort='one', sex='female', y=c(4.548333,
4.591,
3.138,
2.699,
6.622,
6.8795,
5.5925,
1.6715,
4.92775,
6.68525,
4.25775,
8.677))
d3 <- data.frame(cohort='two', sex='male', y=c(7.9645,
16.252,
15.30175,
8.66325,
15.6935,
16.214,
4.056,
8.316,
17.95725,
13.644,
15.76475))
d4 <- data.frame(cohort='two', sex='female', y=c(11.2865,
22.22775,
18.00466667,
12.80925,
16.15425,
14.88133333,
12.0895,
16.5335,
17.68925,
15.00425,
12.149))
d <- rbind(d1, d2, d3, d4)
dd <- datadist(d); options(datadist='dd')

# Fit the default ordinal model (prop. odds)
f <- orm(y ~ cohort * sex, data=d)
f
anova(f)

# Show intercepts as a function of y to estimate the underlying
# conditional distribution.  Result: more uniform than Gaussian
alphas <- coef(f)[1 : num.intercepts(f)]
yunique <- f$yunique[-1]
par(mfrow=c(1,2))
plot(yunique, alphas)

# Compare to distribution of residuals
plot(ecdf(resid(ols(y ~ cohort * sex, data=d))), main='')

M <- Mean(f)
# Confidence intervals for means are approximate
# Confidence intervals for odds ratios or exceedance probabilities
# are correct for ordinal models
Predict(f, cohort, sex, fun=M)

with(d, summarize(y, llist(cohort, sex), smean.cl.normal))
