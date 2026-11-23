require(rms)
require(rmsb)
# require(VGAM) - No - vglm ran out of memory at 24GB with 999 intercepts

# Function to fit rmsb::blrm model that is the same as the orm one, return beta and tau
bfit <- function() {
  ylev <- sort(unique(y))
  # ybar <- mean(ylev[-1])
  g <-  blrm(y ~ x, ppo = ~x, cppo = function(yy) yy - ybar,
             method = "optimizing", iprior = 1, normcppo = FALSE,
    data=data.frame(x, y))
  cof <- coef(g)
  cof <- c(cof['x'], tail(cof, 1))
  names(cof) <- c('beta', 'tau')
  cof
}

n <- 1000
ident <- function(y) y - ybar
betatau <- function(fit) {
  cof <- coef(fit)
  l <- length(cof)
  cof[c(l-1, l)]
}
for(i in 1:10) {
  cat('Seed:', i, '\n\n')
  set.seed(i)
  x <- runif(n); y <- runif(n); ybar <- mean(y)
  f <- orm(y ~ x, ~x, ytrans=ident, trace=1)
  cat('blrm MLE:', bfit(), '\n')
  if(f$fail) cat('orm fit failed\n') else {
    cat('orm beta, tau:', betatau(f), '\n')
  }
}
