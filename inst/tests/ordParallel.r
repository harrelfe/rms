require(rms)
require(ggplot2)
getHdata(nhgh)
set.seed(1)
nhgh$ran <- runif(nrow(nhgh))
f <- orm(gh ~ rcs(age, 4) + ran, data=nhgh, x=TRUE, y=TRUE)
ordParallel(f)
dd <- datadist(nhgh); options(datadist='dd')
ordParallel(f, terms=TRUE)
d <- ordParallel(f, maxcuts=30, onlydata=TRUE)
dd2 <- datadist(d); options(datadist='dd2')  # needed for plotting
g <- orm(Yge_cut ~ (age + ran) * rcs(Ycut, 4), data=d, x=TRUE, y=TRUE)
h <- robcov(g, d$obs)
anova(h)
# Plot inter-quartile-range (on linear predictor "terms") age
# effect vs. cutoff y
qu <- quantile(d$age, c(1, 3)/4)
qu
cuts <- sort(unique(d$Ycut))
cuts
z <- contrast(h, list(age=qu[2], Ycut=cuts),
                 list(age=qu[1], Ycut=cuts))
z <- as.data.frame(z[.q(Ycut, Contrast, Lower, Upper)])
ggplot(z, aes(x=Ycut, y=Contrast)) + geom_line() +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2)

