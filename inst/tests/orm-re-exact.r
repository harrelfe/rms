## Exact (non-approximate) marginal maximum likelihood estimation for
## a proportional-odds (cumulative logit) model with a normally
## distributed random intercept, fit by direct numerical integration
## of each cluster's marginal likelihood -- no Laplace or adaptive
## Gauss-Hermite approximation of the random effect -- followed by
## general-purpose optimization.
##
## Model: for cluster c with random intercept u_c ~ N(0, sigma^2),
##   logit(Pr(Y_ci >= y)) = alpha_y + beta*drug_ci + u_c
## for ordered thresholds alpha_1 > alpha_2 > ... > alpha_k (k = number
## of distinct y values minus 1). This works because every cluster here
## has exactly 2 observations, making the integral over u_c cheap and
## accurate to integrate()'s own tolerance -- for larger cluster sizes
## this approach still works in principle but becomes progressively
## more expensive (the integrand becomes a product over more terms,
## and the numerical integration itself may need more care).
##
## Base R only: integrate() for the exact per-cluster marginal
## likelihood, optim() (BFGS) for maximization. No dependence on rms,
## orm.fit, or any AGQ-based machinery -- useful as an independent
## check on those.

d <- structure(list(drug = c("A", "A", "A", "A", "A", "A", "A", "A",
"A", "A", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B"),
    y = c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0, 2, 1.9,
    0.8, 1.1, 0.1, -0.1, 4.4, 5.5, 1.6, 4.6, 3.4), id = structure(c(1L,
    2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 1L, 2L, 3L, 4L, 5L,
    6L, 7L, 8L, 9L, 10L), levels = c("1", "2", "3", "4", "5",
    "6", "7", "8", "9", "10"), class = "factor")), class = "data.frame",
row.names = c(NA, -20L))

## ---- Set up the ordinal (rank-based) response ----
vals <- sort(unique(d$y))
k    <- length(vals) - 1L        # number of intercepts
n    <- nrow(d)
yy   <- match(d$y, vals) - 1L    # 0-indexed rank of each observation
x    <- as.integer(d$drug == "B")
cl   <- as.integer(d$id)

## ---- Threshold reparametrization: guarantees alpha strictly ----
## descending, so any real-valued par is automatically a valid model
thresholds <- function(par_a) {
  a <- numeric(k)
  a[1] <- par_a[1]
  for(j in 2 : k) a[j] <- a[j - 1] - exp(par_a[j])
  a
}

## ---- Per-observation probability given the random intercept u ----
## Standard cumulative-logit convention:
##   Pr(Y=0) = 1 - plogis(alpha_1 + lp)
##   Pr(Y=k) =     plogis(alpha_k + lp)
##   Pr(Y=y) =     plogis(alpha_y + lp) - plogis(alpha_{y+1} + lp),  0<y<k
obs_prob <- function(yi, lp, a) {
  if(yi == 0) return(1 - plogis(a[1] + lp))
  if(yi == k) return(plogis(a[k] + lp))
  plogis(a[yi] + lp) - plogis(a[yi + 1] + lp)
}

cluster_obs <- split(seq_len(n), cl)

## ---- Negative log-likelihood: exact marginal likelihood via ----
## direct numerical integration over each cluster's own u_c
negloglik <- function(par) {
  a     <- thresholds(par[1 : k])
  beta  <- par[k + 1]
  sigma <- exp(par[k + 2])
  if(any(! is.finite(a)) || sigma <= 0) return(1e10)
  tot <- 0
  for(obs in cluster_obs) {
    lp0 <- beta * x[obs]
    integrand <- function(u) {
      sapply(u, function(uu) {
        prod(mapply(function(yi, l) obs_prob(yi, l + uu, a), yy[obs], lp0)) *
          dnorm(uu, 0, sigma)
      })
    }
    val <- tryCatch(
      integrate(integrand, lower = -Inf, upper = Inf,
                rel.tol = 1e-10, abs.tol = 1e-12)$value,
      error = function(e) 1e-300)
    if(val <= 0 || ! is.finite(val)) val <- 1e-300
    tot <- tot + log(val)
  }
  -tot
}

## ---- Starting values: proportion-based thresholds, beta=0, sigma=1 ----
ncum   <- rev(cumsum(table(factor(yy, levels = 0 : k))))
pp     <- pmin(pmax(ncum[-1] / n, 1e-3), 1 - 1e-3)
alpha0 <- sort(qlogis(pp), decreasing = TRUE)
par_a0 <- numeric(k); par_a0[1] <- alpha0[1]
for(j in 2 : k) par_a0[j] <- log(max(alpha0[j - 1] - alpha0[j], 1e-3))
par0 <- c(par_a0, 0, log(1))

## ---- Optimize ----
fnm <- 'orm-re-exact-opt.rds'
if(file.exists(fnm))
  opt <- readRDS(fnm) else {
  opt <- optim(par0, negloglik, method = "BFGS", hessian = TRUE,
               control = list(maxit = 2000, reltol = 1e-12))
  saveRDS(opt, fnm)
}

if(opt$convergence != 0)
  warning("optim() did not report clean convergence (code ", opt$convergence, ")")

alpha_hat <- thresholds(opt$par[1 : k])
beta_hat  <- opt$par[k + 1]
sigma_hat <- exp(opt$par[k + 2])

## ---- Standard errors (delta method for sigma = exp(log sigma)) ----
vc          <- solve(opt$hessian)
se_beta     <- sqrt(vc[k + 1, k + 1])
se_logsigma <- sqrt(vc[k + 2, k + 2])
se_sigma    <- sigma_hat * se_logsigma

## ---- Finite-difference gradient check at the optimum ----
eps  <- 1e-5
grad <- sapply(seq_along(opt$par), function(j) {
  p1 <- opt$par; p1[j] <- p1[j] + eps
  p2 <- opt$par; p2[j] <- p2[j] - eps
  (negloglik(p1) - negloglik(p2)) / (2 * eps)
})

cat("Exact (numerically integrated) marginal ML estimates\n")
cat("======================================================\n")
cat(sprintf("beta (drug B vs A):          %10.6f   SE: %9.6f\n", beta_hat, se_beta))
cat(sprintf("sigma (random intercept SD): %10.6f   SE: %9.6f\n", sigma_hat, se_sigma))
cat(sprintf("-2 log likelihood: %.6f\n", 2 * opt$value))
cat(sprintf("optim() convergence code: %d (0 = converged)\n", opt$convergence))
cat(sprintf("max|finite-difference gradient| at optimum: %.3e\n", max(abs(grad))))
cat("\nThresholds (alpha):\n")
names(alpha_hat) <- paste0("y>=", vals[-1])
print(alpha_hat)

# Check with orm

require(rms)
f <- orm(y ~ drug + cluster(id), data=d)
print(f, intercepts=TRUE)

require(ordinal)
g <- clmm2(factor(y) ~ drug, data=d, random=id, nAGQ=11, Hess=TRUE)
summary(g)

