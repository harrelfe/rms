set.seed(7)

## true params (unconstrained/eta scale)
kstar_true <- 0.7   # k = plogis(0.7) ~ 0.668
a_true     <- -0.9  # baseline rate exp(-0.9) ~ 0.407 at m=0
b_true     <- 0.06  # rate grows mildly with mean time

k_true <- plogis(kstar_true)
cat("True floor k =", k_true, "\n")

n_subj <- 200
times_list <- replicate(n_subj, sort(sample(seq(0, 12, by = 1), size = sample(4:7,1))), simplify = FALSE)

simOne <- function(tt) {
  n <- length(tt)
  lag  <- abs(outer(tt, tt, "-"))
  mid  <- outer(tt, tt, "+") / 2
  rate <- exp(a_true + b_true*mid)
  R <- k_true + (1-k_true)*exp(-rate*lag); diag(R) <- 1
  R <- .nearestCorPD(R)
  L <- chol(R)
  eps <- as.vector(t(L) %*% rnorm(n))
  data.frame(time = tt, resid = eps)
}

dat_list <- lapply(seq_len(n_subj), function(i) {
  d <- simOne(times_list[[i]]); d$subject <- i; d
})
dat <- do.call(rbind, dat_list)
dat$subject <- factor(dat$subject)
dat$y <- 2 + 0.1*dat$time + dat$resid

cat("N obs:", nrow(dat), "N subjects:", n_subj, "\n\n")

fit <- gls(y ~ time, data = dat,
           correlation = corFloorExp(form = ~ time | subject),
           control = glsControl(msMaxIter = 300, tolerance = 1e-6))

cat("--- Fitted correlation params (k, a, b) ---\n")
print(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
cat("\nTrue: k=", k_true, " a=", a_true, " b=", b_true, "\n")

cat("\n--- Fixed effects ---\n")
print(summary(fit)$tTable)
cat("\nlogLik:", logLik(fit), "\n")

cat("\n=== intervals() ===\n")
print(tryCatch(intervals(fit, which = "var-cov"), error = function(e) conditionMessage(e)))

cat("\n=== combine with lme() (random intercept) ===\n")
fit_lme <- lme(y ~ time, random = ~ 1 | subject, data = dat,
               correlation = corFloorExp(form = ~ time | subject),
               control = lmeControl(msMaxIter = 300, opt = "optim"))
print(coef(fit_lme$modelStruct$corStruct, unconstrained = FALSE))
cat("logLik lme:", logLik(fit_lme), "\n")

cat("\n=== edge case: singleton time per subject ===\n")
dat2 <- rbind(dat, data.frame(time = 3, resid = 0, subject = factor(9999), y = 2.3))
fit2 <- tryCatch(
  gls(y ~ time, data = dat2, correlation = corFloorExp(form = ~ time | subject)),
  error = function(e) conditionMessage(e))
cat(if (is.character(fit2)) paste("ERROR:", fit2) else "OK, singleton subject handled\n")

cat("\n=== boundary check at MLE: corr(0)=1 exactly, corr(large lag) -> k ===\n")
cs <- coef(fit$modelStruct$corStruct, unconstrained = FALSE)
kk <- cs["k"]; aa <- cs["a"]; bb <- cs["b"]
Lgrid <- c(0, 1, 5, 12, 50, 500)
m0 <- 6
rate0 <- exp(aa + bb*m0)
print(kk + (1-kk)*exp(-rate0*Lgrid))

cat("\n=== AIC vs stationary alternatives ===\n")
fit_cs  <- gls(y ~ time, data = dat, correlation = corCompSymm(form = ~ time | subject))
fit_exp <- gls(y ~ time, data = dat, correlation = corExp(form = ~ time | subject))
print(AIC(fit, fit_cs, fit_exp))
