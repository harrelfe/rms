# From Dimitris Rizopoulos 2024-02-20
require(rms)
 
gr <- function (params, x, y) {
  # in 'params' first the alpha then the beta
  p <- ncol(x)
  k <- length(params) - p
  ia <- seq_len(k)
  xb <- as.vector(x %*% params[-ia])
  alpha <- params[ia]
  ealpha <- exp(alpha[-1L])
  alpha <- c(-1e3, cumsum(c(alpha[1L], ealpha)), 1e3)
  P1 <- plogis(alpha[y + 1] + xb)
  P2 <- plogis(alpha[y + 2] + xb)
  Q <- P1 - P2
  pq1Q <- P1 * (1.0 - P1) / Q
  pq2Q <- P2 * (1.0 - P2) / Q
  jacobian <- function (etheta) {
    k <- length(etheta) + 1
    mat <- matrix(0.0, k, k)
    mat[, 1L] <- rep(1.0, k)
    for (i in 2L:k) mat[i:k, i] <- etheta[i - 1]
    mat
  }
  gr_alpha <- rowsum.default(pq1Q, y)[-1L] - rowsum.default(pq2Q, y)[-(k + 1L)]
  - c(crossprod(gr_alpha, jacobian(ealpha)), colSums(x * (pq1Q - pq2Q)))
}
 
lL <- function (params, x, y) {
  p <- ncol(x)
  k <- length(params) - p
  ia <- seq_len(k)
  xb <- as.vector(x %*% params[-ia])
  alpha <- params[ia]
  alpha <- c(-1e3, cumsum(c(alpha[1L], exp(alpha[-1L]))), 1e3)
  -sum(log(plogis(alpha[y + 2] + xb) - plogis(alpha[y + 1] + xb)))
}
 
M <- 100L
time_optim <- time_orm <- numeric(M)
deviance_optim <- deviance_orm <- rep(NA_real_, M)
for (m in seq_len(M)) {
  cat('m=', m, '\n')
  set.seed(2025L + m)
  n <- 1000; p <- 50; k <- n - 1
  x <- matrix(rnorm(n * p), nrow = n)
  y <- order(rnorm(n)) - 1
  alpha <- runif(k, -10, -5)
  alpha_ord <- rev(c(alpha[1], cumsum(exp(alpha[-1]))))
  beta  <- runif(p, -0.5, 0.5)
  #####
  scales <- c(1, rep(0.1, k - 1), rep(1, p))
  time_optim[m] <- system.time({
    res_optim <-
      optim(c(alpha, beta), lL, gr, method = "BFGS", x = x, y = y,
            control = list(reltol = 1e-9, parscale = scales, maxit = 150L))
  })['elapsed']
  deviance_optim[m] <- 2 * res_optim$value
  #####
  time_orm[m] <- system.time({
    # res_orm <-
    #   orm.fit(x, y, initial = c(alpha_ord, beta), eps = 1e-9, compstats = FALSE)
    res_orm <-
      orm.fit(x, y, eps = 1e-9, compstats = FALSE)

  })['elapsed']
  if(res_orm$fail) {
    w <- list(x=x, y=y, init=c(alpha_ord, beta))
    saveRDS(w, '/tmp/w.rds')
    stop()
  }
  if (!is.null(res_orm$deviance)) deviance_orm[m] <- res_orm$deviance[2L]
}
 
# how many times orm() did not converge
na_ind <- is.na(deviance_orm)
sum(na_ind)
 
# deviance differences
summary(deviance_optim[!na_ind] - deviance_orm[!na_ind])
 
# timing differences
summary(time_optim[!na_ind] - time_orm[!na_ind])

if(FALSE) {
  require(rms)
  w <- readRDS('/tmp/w.rds')
  f <- orm.fit(w$x, w$y, eps=1e-9, trace=2)
}
