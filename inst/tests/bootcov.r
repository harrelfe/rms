# From Max Gordon <max@gforge.se>
require(rms)
set.seed(1)
center <- factor(sample(letters[1:8],500,TRUE))
treat  <- factor(sample(c('a','b'),  500,TRUE))
y      <- 8*(treat=='b') + rnorm(500,100,20)
f <- ols(y ~ treat*center, x=TRUE, y=TRUE)
g <- bootcov(f, B=50)

lc <- levels(center)
contrast(f, list(treat='b', center=lc),
         list(treat='a', center=lc))
