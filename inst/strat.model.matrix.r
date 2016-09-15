require(rms)
d <- expand.grid(a=c('a1','a2'), b=c('b1','b2'))
d$y <- Surv(c(1,3,2,4))
f <- y ~ a * strat(b)
m <- model.frame(f, data=d)
Terms <- terms(f, specials='strat', data=d)
specials <- attr(Terms, 'specials')
temp <- survival:::untangle.specials(Terms, 'strat', 1)
Terms <- Terms[- temp$terms]
# X <- rms:::Design(m)
# atr <- attr(X, 'Design')
# atr$colnames

model.matrix(Terms, m)
colnames(model.matrix(Terms, m)[, -1, drop=FALSE])

cph(f, data=d)
