require(rms)
set.seed(1)
age <- rnorm(500, 50, 10)
sex <- sample(c('female', 'male'), 500, TRUE)
y   <- sample(0:4, 500, TRUE)
d   <- expand.grid(age=50, sex=c('female', 'male'))

w <- impactPO(y ~ age + sex, nonpo = ~ sex, newdata=d)
w
# Note that PO model is a better model than multinomial (lower AIC)
# since multinomial model's improvement in fit is low in comparison
# with number of additional parameters estimated.  The same is true
# in comparison to the partial PO model.

# Reverse levels of y so stacked bars have higher y located higher
revo <- function(z) {
  z <- as.factor(z)
  factor(z, levels=rev(levels(as.factor(z))))
}

ggplot(w$estimates, aes(x=method, y=Probability, fill=revo(y))) +
  facet_wrap(~ sex) + geom_col() +
  xlab('') + guides(fill=guide_legend(title=''))

d <- expand.grid(sex=c('female', 'male'), age=c(40, 60))
w <- impactPO(y ~ age + sex, nonpo = ~ sex, newdata=d)
w
ggplot(w$estimates, aes(x=method, y=Probability, fill=revo(y))) +
  facet_grid(age ~ sex) + geom_col() +
  xlab('') + guides(fill=guide_legend(title=''))
