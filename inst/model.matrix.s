d <- data.frame(x1=sample(c('a','b','c'), 20, TRUE),
                x2=sample(c('A','B'),     20, TRUE), y=1:20)
f <- y ~ x1 * strat(x2)
strat <- function(x) x
Terms <- terms(f, specials='strat', data=d)
specials <- attr(Terms, 'specials')
stra <- specials$strat
require(survival) # version 2.37-4; has untangle.specials
z <- untangle.specials(Terms, 'strat', 1)
Terms.ns <- Terms[-z$terms]
X <- model.frame(f, data=d)
colnames(model.matrix(Terms.ns, X))

# If don't remove specials before model.matrix, try to do afterward
x <- model.matrix(Terms, X)
colnames(x)
asg <- attr(x, 'assign')
i <- colnames(x) == '(Intercept)' | (grepl('strat\\(', colnames(x)) & !grepl(':', colnames(x)))
x <- x[, !i]
# How to reconstruct asg?
