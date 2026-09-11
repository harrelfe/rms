## orm.rfit.r -- random-intercept fitting internals for orm.fit().
##
## orm.fit() (in orm_fit.s) is now the sole user-facing entry point
## for both fixed-effects and random-intercept ordinal regression --
## supplying `cluster=` dispatches to ormrfit() (defined in this file)
## after reusing orm.fit()'s own existing initial-value refinement
## machinery (built to ignore clustering) as a warm start. The
## standalone orm.rfit() front-end wrapper that used to live in this
## file has been retired; its data-preparation logic (Ocens handling,
## recode2integer, initial values) was never duplicated here to begin
## with -- it now runs exactly once, inside orm.fit() itself, shared
## by both the clustered and unclustered paths.
##
## Claude Sonnet 5 2026-08-30

## Gauss-Hermite quadrature nodes and weights via the Golub-Welsch
## algorithm (eigen-decomposition of the tridiagonal Jacobi matrix for
## Hermite polynomials). Base R, no dependency. Validated against
## brute-force integrate() during development (see design notes).
## Claude Sonnet 5 2026-08-30          entire function
gauss_hermite_quad <- function(n) {
  if(n == 1L) return(list(nodes=0, weights=sqrt(pi)))
  i <- 1 : (n - 1)
  b <- sqrt(i / 2)
  J <- matrix(0, n, n)
  J[cbind(2 : n, 1 : (n - 1))] <- b
  J[cbind(1 : (n - 1), 2 : n)] <- b
  e <- eigen(J, symmetric=TRUE)
  list(nodes=e$values, weights=sqrt(pi) * e$vectors[1, ]^2)
}

## Per-cluster Newton-Raphson mode-finding for the random intercepts,
## vectorized across all clusters simultaneously via one ormeta call
## per Newton step (not one call per cluster). ia/ia2/sgn/ib/nb are
## from ormidx (data-dependent only, computed once by the caller).
## base_lp must be the FULL linear predictor excluding only the random
## intercept itself, i.e. offset + x %*% beta (ormeta has no x/beta
## arguments -- it consumes a complete linear predictor directly).
## Claude Sonnet 5 2026-09-06          entire function
## wtre added: the per-observation weight w(tre_ij) multiplying the
## random effect's own contribution to eta (eta = ... + w(tre)*u_i).
## Defaults to an all-ones vector, exactly reproducing the unweighted
## formulas below at the cost of one harmless vector multiply per
## Newton iteration -- deliberately not gated behind a separate code
## path, since that multiply is negligible next to the ormeta call
## that already dominates this loop's cost.
## Chain rule: d(eta_ij)/d(u_i) = w(tre_ij), so the per-cluster score
## contribution scales by w (not 1) and the curvature by w^2 (not 1),
## exactly as for any linear reparametrization of a scalar argument.
## Claude Sonnet 5 2026-09-08          entire function
## Extended to implement step-halving on ormeta's salloc=999, per its
## own documented contract ("signals the caller to step-halve, same
## as ormll's what=1 path") -- confirmed directly, from ormeta.f90's
## actual source, that this was NOT previously honored here: any
## nonzero salloc was treated as an immediate, unrecoverable failure,
## even though it is a designed "please step-halve" signal, not a
## hard error. This was found to matter in practice specifically
## under the probit link: its CDF (via erf()) saturates to exactly
## 0/1 at a much smaller |eta| than the logistic CDF does, so a
## Newton step that overshoots into a region where two adjacent
## intercepts' cumulative probabilities both round to the same
## saturated value (making their difference, a category probability,
## exactly zero) is considerably more likely to occur than under
## logistic -- where the ordinary (non-clustered) Newton-Raphson loop
## already step-halves through the identical condition via ormll,
## which is why a probit fit without `cluster` can succeed cleanly
## while the random-intercept extension, missing this same handling,
## failed outright.
clusterModeFind <- function(alpha, beta, sigma, gamma_c_init, base_lp,
                            wt, ia, ia2, sgn, ib, nb, cluster, nc,
                            link, k, n, wtre=rep(1e0, n), maxit=30L, tol=1e-8,
                            minstepsize=1e-4) {
  ## Claude Sonnet 5 2026-09-08
  ## deb() is a real print-with-label function (see Hmisc::Fdebug) only
  ## when options(orm.re.debug=TRUE) is set; otherwise it's a do-nothing
  ## function costing almost no time (its argument, however expensive,
  ## is never forced) -- replaces the earlier ad-hoc "always retry the
  ## one specific failing call with debug=1" logic with a general-
  ## purpose mechanism usable at any point in this function, not just
  ## that one call.
  deb <- Hmisc::Fdebug('orm.re.debug')
  gamma_c <- gamma_c_init
  ## Guard against non-finite lp BEFORE calling .Fortran: R's own
  ## .Fortran() interface performs a hard, unconditional check for
  ## NA/NaN/Inf in numeric arguments and throws an R-level error if it
  ## finds one, BEFORE the Fortran code (and hence salloc) is ever
  ## reached -- confirmed directly to happen in practice (a large,
  ## still-growing sigma can drive curv toward its -1e-8 floor,
  ## producing an enormous, eventually non-finite, Newton step).
  ## Treated as an "invalid trial point" exactly like salloc=999 --
  ## returning a synthetic salloc lets the SAME step-halving loop
  ## below handle both cases uniformly.
  eval_ormeta <- function(gc) {
    gamma_obs <- wtre * gc[cluster]
    lp <- base_lp + gamma_obs
    deb(alpha); deb(gc); deb(lp)
    if(any(! is.finite(lp))) return(list(salloc=999L))
    w <- .Fortran(F_ormeta, n, k, link, alpha, lp1=lp, lp2=lp, wt,
                 ia, ia2, sgn, ib, nb,
                 logd=numeric(n), g=numeric(n), h=numeric(n),
                 s1=numeric(n), s2=numeric(n),
                 debug=0L, salloc=integer(1))
    deb(w$salloc)
    w
  }
  ## No prior, known-good point exists to step-halve back toward here
  ## -- this can only fail if (alpha, beta) ALONE, with no random
  ## effect contribution yet (gamma_c_init is 0 at the very start of
  ## fitting), already produces an exactly-zero category probability
  ## for some observation -- a different, more fundamental problem
  ## step-halving cannot address. Confirmed to occur under probit even
  ## when the SAME (alpha, beta) converged cleanly via ormll in the
  ## ordinary (non-clustered) fit -- plausibly a category probability
  ## that is positive but vanishingly small there, tipped to exactly
  ## zero here by a floating-point-level difference between R's own
  ## x%*%beta and ormll's internal computation of the same quantity,
  ## a difference probit's much-faster-saturating CDF makes far more
  ## consequential than logistic's would. Set options(orm.re.debug=TRUE)
  ## to see every gc/lp/salloc value printed above, when this or any
  ## other failure here needs diagnosing.
  w <- eval_ormeta(gamma_c)
  if(w$salloc != 0) return(list(fail=TRUE, code=w$salloc))

  for(it in 1 : maxit) {
    Sg    <- as.vector(rowsum(wtre * w$g, cluster, reorder=TRUE))
    Sh    <- as.vector(rowsum(wtre^2 * w$h, cluster, reorder=TRUE))
    grad  <- Sg - gamma_c / sigma^2
    curv  <- Sh - 1 / sigma^2        # negative at a proper maximum
    curv  <- pmin(curv, -1e-8)       # guard against numerical breakdown
    step  <- grad / curv

    step_size <- 1.0
    repeat {
      gamma_c_new <- gamma_c - step_size * step
      w_new <- eval_ormeta(gamma_c_new)
      if(w_new$salloc == 0) break
      step_size <- step_size / 2
      if(step_size < minstepsize) return(list(fail=TRUE, code=w_new$salloc))
    }
    w <- w_new

    if(max(abs(gamma_c_new - gamma_c)) < tol) { gamma_c <- gamma_c_new; break }
    gamma_c <- gamma_c_new
  }
  list(fail=FALSE, gamma_c=gamma_c, kappa_c=-curv, iter=it)
}

## AGQ log-likelihood (and, if score=TRUE, the accumulated score/Hessian
## for alpha, beta, & sigma1/sigma2) at given (alpha, beta, sigma), using a
## pre-computed mode (gamma_c, kappa_c). ormeta supplies the cheap
## per-node cluster log-likelihoods; ormll (called once per node with
## reweighted wt and node-shifted offset) supplies the score/Hessian,
## reusing its existing, verified accumulation machinery unchanged.
## The penalty matrix is passed as zero to each per-node ormll call
## and applied exactly once afterward, to avoid counting it nAGQ times.
##
## base_lp here means offset + x %*% beta (full, as in clusterModeFind)
## when building ormeta's lp1/lp2; when building ormll's own `offset`
## argument, x %*% beta must be subtracted back out first, since ormll
## adds it itself internally -- see node_off below.
##
## Claude Sonnet 5 2026-09-06          entire function
## Extended for the mre (multiplier for random effects) weighting
## feature: eta = ... + [sigma1*(1-mre) + sigma2*mre]*v_i, where v_i
## ~ N(0,1) is a STANDARDIZED per-cluster random effect (unlike the
## plain, mre-inactive case below, where gamma_c directly represents
## u_i ~ N(0,sigma^2) with sigma itself estimated). mre is a FIXED,
## user-specified quantity (0 by default at each subject's anchor
## observation, 1 thereafter, for the recommended Markov-1 usage);
## sigma1 and sigma2 are the two estimated scale parameters, entering
## ADDITIVELY (each multiplying its own fixed, known pseudo-covariate)
## rather than one multiplicatively modifying the other.
##
## This design replaces two earlier attempts, both found wanting on
## direct testing:
## - (wl, rho): w(tre)=wl+(1-wl)*rho^tre entered NONLINEARLY and was
##   genuinely unstable -- an early, poorly-conditioned step could
##   drive rho to its exact floating-point boundary, collapsing wl's
##   curvature; a direct grid-scan then showed this was a genuine
##   competing-mode/weak-identifiability property of that
##   parametrization, not a fixable numerical issue.
## - (sigma, re_mult): (1-mre*re_mult)*u_i, u_i~N(0,sigma^2), entered
##   LINEARLY and was well-identified in isolation, but re_mult's
##   scale is defined only RELATIVE to sigma -- as sigma->0, re_mult
##   becomes structurally meaningless and can drift without bound,
##   confirmed directly (a true re_mult>1 scenario diverged: sigma
##   collapsed toward 0 while re_mult diverged toward -infinity, even
##   though the true parameter values gave a clearly higher
##   likelihood than where the fit was heading).
## sigma1/sigma2 avoid this: each is an ABSOLUTE scale, identified
## directly by the variance actually present in its own share of the
## data (weighted by (1-mre) or mre respectively), with no ratio
## relationship between them. Verified directly via curvature checks
## -- including deliberately small samples (nc=22), a true sigma2
## implying a NEGATIVE net weight, and genuinely continuous (not just
## two-point) mre values spread outside [0,1] -- correlation between
## sigma1 and sigma2 stayed under 0.12 in magnitude throughout, a
## dramatic improvement over (sigma, re_mult)'s 0.7-0.8.
##
## sigma1 is constrained positive (the joint Newton step estimates
## log(sigma1)) to pin down v_i's otherwise-arbitrary sign convention
## (flipping (sigma1,sigma2,v_i) -> (sigma1,-sigma2,-v_i) leaves the
## distribution unchanged unless something fixes the sign); sigma2 is
## left unconstrained to preserve the "subtract/reverse" flexibility
## that motivated this whole feature.
##
## mre must have within-cluster variation for sigma1/sigma2 to be
## separately identifiable at all -- see ormrfit's own validation
## (diff(range(mre))==0 is rejected) for the degenerate cases this
## rules out (mre constant at 0, at 1, or at any other single value).
agqStep <- function(alpha, beta, sigma, gamma_c, kappa_c, nAGQ,
                    base_lp, x, wt, ia, ia2, sgn, ib, nb,
                    cluster, nc, link, k, p, n, penmat, penhess=1L,
                    intcens=0L, score=TRUE, mre=NULL, sigma1=NULL, sigma2=NULL) {
  gh     <- gauss_hermite_quad(nAGQ)
  step_c <- sqrt(2 / kappa_c)

  use_mre <- length(mre) > 0
  ## w = sigma1*(1-mre) + sigma2*mre: fixed per observation given the
  ## current sigma1/sigma2, the SAME at every AGQ node (only
  ## gamma_node_c itself varies by node, not its weight) -- computed
  ## once, outside the node loop. NOTE: when use_mre, the CALLER is
  ## responsible for passing sigma=1 (v_i's fixed, standardized prior
  ## SD) -- sigma1/sigma2, not this function's own `sigma` argument,
  ## carry the actual estimated scale information in that case.
  wtre <- if(use_mre) sigma1 * (1 - mre) + sigma2 * mre else rep(1e0, n)

  lt <- matrix(0e0, nc, nAGQ)
  for(m in 1 : nAGQ) {
    gamma_node_c <- gamma_c + step_c * gh$nodes[m]
    gamma_node   <- wtre * gamma_node_c[cluster]
    lp       <- base_lp + gamma_node
    ## Claude Sonnet 5 2026-09-08          2 lines
    ## Same guard as clusterModeFind: R's .Fortran() throws an
    ## uncatchable hard error on non-finite arguments before any
    ## salloc-based handling could run -- confirmed directly to be
    ## reachable here too (e.g. a large, still-growing sigma widening
    ## step_c enough to put some AGQ node's gamma_node_c at an extreme
    ## value).
    if(any(! is.finite(lp))) return(list(fail=TRUE, code=999L))
    w        <- .Fortran(F_ormeta, n, k, link, alpha, lp1=lp, lp2=lp, wt,
                         ia, ia2, sgn, ib, nb,
                         logd=numeric(n), g=numeric(n), h=numeric(n),
                         s1=numeric(n), s2=numeric(n),
                         debug=0L, salloc=integer(1))
    if(w$salloc != 0) return(list(fail=TRUE, code=w$salloc))
    ## Claude Sonnet 5 2026-08-30          1 line
    ## The Gaussian random-effect prior log-density is part of the
    ## integrand being quadrature-approximated and must be added here
    ## (clusterModeFind already includes it in the mode/curvature via
    ## grad = Sg - gamma_c/sigma^2 etc.; omitting it here would leave sigma
    ## completely unpenalized and the "likelihood" unbounded in sigma).
    lt[, m] <- log(gh$weights[m]) + gh$nodes[m]^2 +
               as.vector(rowsum(w$logd, cluster, reorder=TRUE)) +
               dnorm(gamma_node_c, mean=0, sd=sigma, log=TRUE)
  }

  mrow      <- apply(lt, 1, max)
  rho_w     <- exp(lt - mrow)
  denom     <- rowSums(rho_w)
  loglik_c  <- log(step_c) + mrow + log(denom)   # natural-log scale
  rho_w     <- rho_w / denom                     # nc x nAGQ, rows sum to 1
  total_ll  <- sum(loglik_c)

  if(! score)
    return(list(fail=FALSE, loglik=total_ll, rho=rho_w, gh=gh, step_c=step_c,
               gamma_c=gamma_c, kappa_c=kappa_c))

  xbeta <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)

  nai <- if(intcens == 1L) 1000000L else 1L

  ## p_star: p, plus two more columns for sigma1/sigma2 when mre is
  ## active (replacing the single re_mult column of the earlier design).
  p_star  <- p + if(use_mre) 2L else 0L
  zeromat_star <- matrix(0e0, p_star, p_star)

  accum_grad <- numeric(k + p_star)
  accum_ha   <- if(intcens == 0L) matrix(0e0, k, 2) else NULL
  accum_row  <- NULL; accum_col <- NULL; accum_ai <- NULL; accum_ne <- 0L
  accum_hb   <- matrix(0e0, p_star, p_star)
  accum_hab  <- matrix(0e0, k, p_star)
  ## Claude Sonnet 5 2026-09-06          3 lines
  ## tau=log(sigma)'s own gradient/curvature: only meaningful (and
  ## only accumulated) when mre is INACTIVE -- when mre is active,
  ## sigma is fixed at 1 (v_i's own prior), not estimated at all;
  ## sigma1/sigma2 are handled via the pseudo-covariate columns below
  ## instead, exactly like beta.
  accum_tau_grad <- 0; accum_tau_curv <- 0

  for(m in 1 : nAGQ) {
    gamma_node_c  <- gamma_c + step_c * gh$nodes[m]
    gamma_node    <- wtre * gamma_node_c[cluster]
    node_off  <- base_lp - xbeta + gamma_node        # ormll adds x%*%beta itself
    node_wt   <- wt * rho_w[cluster, m]
    a_arg     <- if(intcens == 0L) matrix(0e0, k, 2) else matrix(0e0, 0, 2)
    ## Claude Sonnet 5 2026-09-08          2 lines
    ## Same non-finite-argument guard as the log-likelihood loop above
    ## (and clusterModeFind) -- see their comments for why this is
    ## needed before .Fortran, here for the score/Hessian call instead.
    if(any(! is.finite(node_off))) return(list(fail=TRUE, code=999L))

    if(! use_mre) {
      ## Claude Sonnet 5 2026-09-06          3 lines
      ## tau's own per-node contribution: d/dtau[log dnorm(u;0,sigma)] =
      ## -1+u^2/sigma^2 (matching sparseMissingInfo's identical "tau
      ## column" formula), second derivative -2*u^2/sigma^2. AGQ-weighted
      ## by this node's posterior cluster weights, exactly like the data
      ## log-lik terms are weighted via node_wt above.
      accum_tau_grad <- accum_tau_grad + sum(rho_w[, m] * (-1 + gamma_node_c^2 / sigma^2))
      accum_tau_curv <- accum_tau_curv + sum(rho_w[, m] * (-2 * gamma_node_c^2 / sigma^2))
    }

    ## Claude Sonnet 5 2026-09-06          12 lines
    ## Pseudo-covariate columns for (log(sigma1), sigma2): eta's random-
    ## effect term is [sigma1*(1-mre)+sigma2*mre]*gamma_node_c[cluster].
    ## d(eta)/d(sigma1) = (1-mre)*gamma_node_c[cluster]; the joint
    ## Newton step estimates log(sigma1) (to keep it positive), so an
    ## extra chain-rule factor of sigma1 (=d(sigma1)/d(log sigma1))
    ## applies: d(eta)/d(log sigma1) = (1-mre)*gamma_node_c[cluster]*
    ## sigma1. d(eta)/d(sigma2) = mre*gamma_node_c[cluster] directly,
    ## no transform (sigma2 is left unconstrained). Passing beta=0 for
    ## both columns is correct, not a placeholder -- same construction
    ## as the (wl,rho) and re_mult pseudo-columns this replaced: eta
    ## already has the correctly-weighted contribution baked into
    ## gamma_node above, so a zero coefficient here leaves eta
    ## unchanged while ormll's own score machinery still returns
    ## exactly the score/Hessian contribution evaluated at that point.
    if(use_mre) {
      d_logsigma1 <- (1 - mre) * gamma_node_c[cluster] * sigma1
      d_sigma2    <- mre * gamma_node_c[cluster]
      x_star    <- cbind(x, d_logsigma1, d_sigma2)
      beta_star <- c(beta, 0, 0)
    } else { x_star <- x; beta_star <- beta }

    wc <- .Fortran(F_ormll, n, k, p_star, x_star, ia, ia2, sgn, ib, nb,
                   node_off, node_wt, zeromat_star,
                   link=link, alpha, beta_star, logL=numeric(1),
                   grad=numeric(k + p_star), lpe=numeric(n),
                   a=a_arg, b=matrix(0e0, p_star, p_star), ab=matrix(0e0, k, p_star),
                   intcens=intcens, row=integer(nai), col=integer(nai), ai=numeric(nai),
                   nai=nai, ne=integer(1),
                   ## nu=0L intentionally: this disables ormll's sparse
                   ## score-matrix computation entirely (nu is both a
                   ## buffer size AND a behavioral switch via ormll's
                   ## own `if(nu > 0)` gate) -- urow/ucol/um are
                   ## unused placeholders in that case, so a length-1
                   ## buffer (larger than the nu=0 the Fortran code
                   ## actually honors) is safe, unlike g/h above which
                   ## are unconditionally written regardless of size.
                   urow=integer(1), ucol=integer(1), um=numeric(1), nu=0L, nuu=integer(1),
                   what=3L, debug=0L, penhess=0L, salloc=integer(1))
    if(wc$salloc != 0) return(list(fail=TRUE, code=wc$salloc))
    accum_grad <- accum_grad + wc$grad
    accum_hb   <- accum_hb   + wc$b
    accum_hab  <- accum_hab  + wc$ab

    if(intcens == 0L) {
      accum_ha <- accum_ha + wc$a
    } else {
      ## Claude Sonnet 5 2026-08-30          8 lines
      ## Verified empirically (design notes, "Step 0") that ormll's
      ## sparse (row, col, ne) index structure for interval-censored
      ## data depends only on ia/ia2 (hence on y/y2/k), never on
      ## alpha/beta/offset/wt -- so it is identical across every AGQ
      ## node, and the sparse Hessian values (ai) can be accumulated
      ## elementwise using the FIRST node's index structure, exactly
      ## like the dense ha case above. Re-checked (not just assumed)
      ## on every subsequent node in case a future change to ormll
      ## ever breaks that invariant.
      if(m == 1) {
        accum_ne  <- wc$ne
        accum_row <- wc$row[1 : accum_ne]
        accum_col <- wc$col[1 : accum_ne]
        accum_ai  <- wc$ai[1 : accum_ne]
      } else {
        if(wc$ne != accum_ne ||
          ! identical(wc$row[1 : wc$ne], accum_row) ||
          ! identical(wc$col[1 : wc$ne], accum_col))
          stop('ormll returned a different sparse intercept-Hessian ',
              'index structure across AGQ nodes -- the assumption this ',
              'accumulation relies on has been violated; see the design ',
              'notes ("Step 0") for orm.rfit interval-censoring support.')
        accum_ai <- accum_ai + wc$ai[1 : accum_ne]
      }
    }
  }

  ## Claude Sonnet 5 2026-09-06          14 lines
  ## Exact second-derivative correction for the log(sigma1)
  ## reparametrization. ormll's own Hessian entry for the
  ## d_logsigma1=(1-mre)*gamma_node_c*sigma1 pseudo-covariate treats it
  ## as a FIXED quantity, capturing only d2L/dtheta2*(dtheta/dphi)^2
  ## (phi=log(sigma1), theta=sigma1=exp(phi)) -- missing the
  ## dL/dtheta*d2(theta)/dphi2 term a TRUE second derivative through
  ## this nonlinear link requires. For theta=exp(phi), d(theta)/dphi =
  ## d2(theta)/dphi2 = theta, so this missing term equals exactly
  ## log(sigma1)'s own (already correctly chain-ruled) gradient.
  ## Confirmed directly against finite differences of the true
  ## surrogate objective (exact agreement) and confirmed as the root
  ## cause of an observed collapse: without it, ormll's Hessian entry
  ## vanishes like sigma1^2 near sigma1=0 (faster than the gradient's
  ## sigma1), so a Newton step (~grad/hess ~ 1/sigma1) grows without
  ## bound and drives sigma1 toward the very boundary that ought to be
  ## resisted -- a self-reinforcing artifact of the approximation, not
  ## a real feature of the likelihood.
  if(use_mre) {
    logsigma1_idx <- p + 1   # position within the p_star-sized beta-like block
    accum_hb[logsigma1_idx, logsigma1_idx] <- accum_hb[logsigma1_idx, logsigma1_idx] +
                                              accum_grad[k + logsigma1_idx]
  }

  ## Penalty applies only to the original p beta columns, never to the
  ## sigma1/sigma2 pseudo-columns -- indexed explicitly since accum_hb
  ## may be larger than p x p.
  if(p > 0) {
    accum_grad[(k + 1) : (k + p)] <- accum_grad[(k + 1) : (k + p)] -
                                     as.vector(penmat %*% beta)
    if(penhess > 0) accum_hb[1 : p, 1 : p] <- accum_hb[1 : p, 1 : p] - penmat
  }

  list(fail=FALSE, loglik=total_ll, rho=rho_w, gh=gh, step_c=step_c,
       gamma_c=gamma_c, kappa_c=kappa_c,
       grad=accum_grad, ha=accum_ha, hb=accum_hb, hab=accum_hab,
       row=accum_row, col=accum_col, ai=accum_ai, ne=accum_ne,
       tau_grad=accum_tau_grad, tau_curv=accum_tau_curv,
       p_star=p_star)
}

## Claude Sonnet 5 2026-09-07
## Re-added (was removed as dead code when sigma was folded into the
## joint Newton step for all cases -- that fold was then found to
## regress the mre-inactive case, see onestep's own revision note).
## Used ONLY when mre is inactive: a separate, direct 1-D profile
## search over sigma, exactly as in the original, pre-fold design.
sigmaObjective <- function(log_sigma, alpha, beta, gamma_c_init, base_lp,
                           wt, ia, ia2, sgn, ib, nb, cluster, nc,
                           link, k, n, nAGQ, maxit.mode) {
  sigma <- exp(log_sigma)
  mf <- clusterModeFind(alpha, beta, sigma, gamma_c_init, base_lp,
                        wt, ia, ia2, sgn, ib, nb, cluster, nc,
                        link, k, n, maxit=maxit.mode)
  if(mf$fail) return(-Inf)
  ag <- agqStep(alpha, beta, sigma, mf$gamma_c, mf$kappa_c, nAGQ,
               base_lp, x=0, wt, ia, ia2, sgn, ib, nb, cluster, nc,
               link, k, p=0L, n, penmat=matrix(0, 0, 0), score=FALSE)
  if(ag$fail) return(-Inf)
  ag$loglik
}

## Internal random-intercept fitting workhorse, called from orm.fit()
## (orm_fit.s) when `cluster` is supplied. Not exported, not meant to
## be called directly -- see orm.fit()'s own documentation.
## Claude Sonnet 5 2026-08-30
## y2      : upper endpoint of interval-censored outcome coding, as in
##           orm.fit()'s internal ormfit; y2=y for uncensored observations
## cluster : integer vector with values 1, 2, ..., nc
## nc      : number of distinct clusters, max(cluster)
##
## Claude Sonnet 5 2026-09-06
## Extended for mre (multiplier for random effects) weighting: the
## per-cluster random effect v_i ~ N(0,1) (standardized -- see
## agqStep's header comment) is scaled by [sigma1*(1-mre)+sigma2*mre]
## at each observation, where mre is a FIXED, user-specified quantity
## and sigma1/sigma2 are estimated. NOTE the shape difference from the
## mre-inactive case below: when mre is active there is no separate
## "tau" (log-sigma) correction at all -- sigma1/sigma2 are ordinary
## score columns (like beta's), since the random effect's own prior SD
## is fixed at 1, not estimated. mre must vary WITHIN at least some
## clusters for sigma1/sigma2 to be separately identifiable -- ormrfit
## validates this directly (diff(range(mre))==0 is rejected, ruling
## out mre constant at 0, at 1, or at any other single value).
sparseMissingInfo <- function(alpha, beta, sigma, gamma_c, kappa_c, base_lp,
                              x, wt, ia, ia2, sgn, ib, nb,
                              cluster, nc, link, k, p, n, nAGQ,
                              mre=NULL, sigma1=NULL, sigma2=NULL) {
  gh     <- gauss_hermite_quad(nAGQ)
  step_c <- sqrt(2 / kappa_c)
  use_mre <- length(mre) > 0
  ## Claude Sonnet 5 2026-09-06          4 lines
  ## p_wt: total extra parameters beyond alpha/beta. mre-active: 2
  ## (log_sigma1, sigma2), no separate tau at all. mre-inactive: 1
  ## (tau=log(sigma), the ONLY case that still needs the tau-tau
  ## complete_tau_info correction below).
  p_wt <- if(use_mre) 2L else 1L
  ptot <- p + p_wt   # p fixed-effect columns plus p_wt random-effect-scale columns
  wtre <- if(use_mre) sigma1 * (1 - mre) + sigma2 * mre else rep(1e0, n)

  obs_by_cluster <- split(seq_len(n), cluster)

  logd_nodes <- matrix(0e0, n, nAGQ)
  s1_nodes   <- matrix(0e0, n, nAGQ)
  s2_nodes   <- matrix(0e0, n, nAGQ)
  g_nodes    <- matrix(0e0, n, nAGQ)
  gamma_node_mat_full <- matrix(0e0, nc, nAGQ)   # per-cluster node values, all nodes

  for(m in 1 : nAGQ) {
    gamma_node_c <- gamma_c + step_c * gh$nodes[m]
    gamma_node_mat_full[, m] <- gamma_node_c
    gamma_node   <- wtre * gamma_node_c[cluster]
    lp       <- base_lp + gamma_node
    ## Claude Sonnet 5 2026-09-08          1 line
    ## Same non-finite-argument guard as clusterModeFind/agqStep --
    ## see their own comments for why this is needed before .Fortran.
    if(any(! is.finite(lp)))
      stop('non-finite linear predictor encountered while computing the ',
          'sparse missing-information correction (sigma or a random ',
          'effect may not have converged to a finite value)')
    w <- .Fortran(F_ormeta, n, k, link, alpha, lp1=lp, lp2=lp, wt,
                 ia, ia2, sgn, ib, nb,
                 logd=numeric(n), g=numeric(n), h=numeric(n),
                 s1=numeric(n), s2=numeric(n),
                 debug=0L, salloc=integer(1))
    if(w$salloc != 0)
      stop('ormeta failed while computing the sparse missing-information correction')
    logd_nodes[, m] <- w$logd
    s1_nodes[, m]   <- w$s1
    s2_nodes[, m]   <- w$s2
    g_nodes[, m]    <- w$g
  }

  ## Per-cluster, per-node AGQ log-lik and posterior node weights.
  ## sd=1 in the dnorm call below when use_mre (v_i's fixed prior);
  ## the passed-in `sigma` argument is only meaningful when !use_mre.
  clusterloglik <- matrix(0e0, nc, nAGQ)
  for(j in 1 : nc) {
    obs <- obs_by_cluster[[j]]
    clusterloglik[j, ] <- colSums(logd_nodes[obs, , drop=FALSE])
  }
  prior_sd <- if(use_mre) 1 else sigma
  lt   <- sweep(clusterloglik, 2, log(gh$weights) + gh$nodes^2, "+") +
          dnorm(gamma_node_mat_full, mean=0, sd=prior_sd, log=TRUE)
  mrow <- apply(lt, 1, max)
  rho_w <- exp(lt - mrow); rho_w <- rho_w / rowSums(rho_w)   # nc x nAGQ

  ## Claude Sonnet 5 2026-09-06          5 lines
  ## Complete-data tau-tau diagonal correction: ONLY meaningful when
  ## mre is inactive (sigma is estimated via a separate tau, not via
  ## ordinary score columns). When mre is active this is unused --
  ## sigma1/sigma2 get their own score columns below, exactly like
  ## beta, with no separate variance-parameter correction needed.
  complete_tau_info <- if(use_mre) 0 else
    sum(rho_w * gamma_node_mat_full^2) * 2 / sigma^2

  ## Per-cluster small dense outer products, accumulated into a
  ## global sparse triplet list. Dimension na + ptot: when mre is
  ## inactive, ptot's last position is tau (log-sigma), whose own
  ## per-node score contribution is -1 + gamma_node^2/sigma^2 -- the
  ## constant -1 contributes nothing to the AGQ-weighted covariance
  ## (already subtracted via Sbarj below) but is included for clarity/
  ## direct correspondence to the math. When mre is active, the last
  ## two positions are log_sigma1/sigma2 instead, treated as ordinary
  ## score columns with no such special-cased entry.
  tI <- integer(0); tJ <- integer(0); tX <- numeric(0)

  for(j in 1 : nc) {
    obs    <- obs_by_cluster[[j]]
    iao    <- ia[obs]; ia2o <- ia2[obs]
    active <- sort(unique(c(iao, ia2o[ia2o > 0])))
    na     <- length(active)
    dimj   <- na + ptot

    Sj <- matrix(0e0, nAGQ, dimj)
    for(m in 1 : nAGQ) {
      if(na > 0) {
        s1o <- s1_nodes[obs, m]; s2o <- s2_nodes[obs, m]
        for(l in 1 : na)
          Sj[m, l] <- sum(s1o[iao == active[l]]) + sum(s2o[ia2o == active[l]])
      }
      if(p > 0) {
        go <- g_nodes[obs, m]
        for(l in 1 : p) Sj[m, na + l] <- sum(go * x[obs, l])
      }
      if(use_mre) {
        ## Claude Sonnet 5 2026-09-06          6 lines
        ## log_sigma1/sigma2 score columns, same pseudo-covariates as
        ## agqStep's d_logsigma1/d_sigma2 -- gamma_node_mat_full[j,m]
        ## is this cluster's own node value.
        go <- g_nodes[obs, m]
        gn <- gamma_node_mat_full[j, m]
        d_logsigma1 <- (1 - mre[obs]) * gn * sigma1
        d_sigma2    <- mre[obs] * gn
        Sj[m, na + p + 1] <- sum(go * d_logsigma1)
        Sj[m, na + p + 2] <- sum(go * d_sigma2)
      } else {
        Sj[m, dimj] <- -1 + gamma_node_mat_full[j, m]^2 / sigma^2   # tau column
      }
    }
    wj    <- rho_w[j, ]
    Sbarj <- as.vector(wj %*% Sj)
    Dj    <- sweep(Sj, 2, Sbarj)
    Vj    <- crossprod(Dj, wj * Dj)   # AGQ-weighted covariance, dimj x dimj

    idxglobal <- c(active,
                  if(p > 0) (k + 1) : (k + p) else integer(0),
                  (k + p + 1) : (k + p + p_wt))
    tI <- c(tI, rep(idxglobal, times = dimj))
    tJ <- c(tJ, rep(idxglobal, each  = dimj))
    tX <- c(tX, as.vector(Vj))
  }

  ptotall <- k + p + p_wt
  missing_info <- if(length(tI) > 0)
    Matrix::sparseMatrix(i=tI, j=tJ, x=tX, dims=c(ptotall, ptotall)) else
    Matrix::Matrix(0, ptotall, ptotall, sparse=TRUE)

  list(missing_info = missing_info, complete_tau_info = complete_tau_info)
}

sparseInfoMatrix <- function(alpha, beta, sigma, gamma_c, kappa_c, base_lp,
                             x, wt, ia, ia2, sgn, ib, nb, cluster, nc,
                             link, k, p, n, ha, hb, hab, row, col, ai, ne,
                             intcens, nAGQ, iname, xname,
                             mre=NULL, sigma1=NULL, sigma2=NULL) {
  ## Claude Sonnet 5 2026-08-30
  ## Branches on intcens for how the complete-data alpha-alpha block
  ## is built: the tri-band ha (k x 2) when intcens=0, or directly
  ## from agqStep's own accumulated sparse (row, col, ai) triplet
  ## when intcens=1 -- the same general triplet form this function's
  ## own output already uses for info$a regardless of intcens (see
  ## the earlier reshaping to match ormfit's format), so no further
  ## reshaping is needed here beyond negating for the info-matrix
  ## sign convention.
  Ha_sparse <- if(intcens == 1L)
    Matrix::sparseMatrix(row, col, x=ai, dims=c(k, k), symmetric=TRUE)
    else if(k > 1)
    Matrix::bandSparse(k, k=c(0, 1),
                       diagonals=list(ha[, 1], ha[1 : (k - 1), 2]),
                       symmetric=TRUE) else
    Matrix::Matrix(ha[1, 1], 1, 1, sparse=TRUE)

  ## Claude Sonnet 5 2026-09-06
  ## hb/hab as passed in are already p_star x p_star / k x p_star,
  ## where p_star=p (mre inactive) or p+2 (mre active, sigma1/sigma2
  ## already included -- jointly Newton-updated with alpha/beta, no
  ## separate profile search or include_wt_cols distinction needed).
  ## When mre is inactive, log(sigma) (tau) is genuinely new on top of
  ## that, extending p_star by exactly one more row/column -- when
  ## mre is active there is NO such extra row/column: sigma1/sigma2
  ## already ARE p_star's extra columns, with no separate tau at all
  ## (the random effect's own prior SD is fixed at 1 in that case).
  ## Complete-data alpha/beta-tau cross terms are exactly zero for the
  ## same reason beta-tau's are (the data likelihood given u doesn't
  ## involve sigma) -- so, when applicable, that extra row/column is
  ## zero except the tau-tau diagonal entry, which sparseMissingInfo
  ## computed alongside everything else it already needed
  ## (complete_tau_info, on the raw Hessian scale, so negated below
  ## along with everything else for the info-matrix sign convention;
  ## complete_tau_info is 0, and unused, when mre is active).
  use_mre <- length(mre) > 0
  p_star  <- p + if(use_mre) 2L else 0L    # matches ag.final$hb/hab's actual dimension
  ptau    <- if(use_mre) p_star else p_star + 1L
  hb_ext    <- matrix(0e0, ptau, ptau)
  if(p_star > 0) hb_ext[1 : p_star, 1 : p_star] <- hb
  hab_ext   <- matrix(0e0, k, ptau)
  if(p_star > 0) hab_ext[, 1 : p_star] <- hab

  mi <- sparseMissingInfo(alpha, beta, sigma, gamma_c, kappa_c, base_lp,
                          x, wt, ia, ia2, sgn, ib, nb,
                          cluster, nc, link, k, p, n, nAGQ,
                          mre=mre, sigma1=sigma1, sigma2=sigma2)
  if(! use_mre) hb_ext[ptau, ptau] <- -mi$complete_tau_info   # raw Hessian scale, negated below with the rest

  complete_data_info <-
    rbind(cbind(-Ha_sparse,                -Matrix::Matrix(hab_ext, sparse=TRUE)),
          cbind(-Matrix::Matrix(t(hab_ext), sparse=TRUE), -Matrix::Matrix(hb_ext, sparse=TRUE)))

  V <- complete_data_info - mi$missing_info

  ## Split into the a (sparse triplet, upper-triangle only) / b / ab
  ## pieces ormfit's own intcens=1 output uses. b is now
  ## ptau x ptau and ab is k x ptau, with the trailing columns being
  ## log(sigma) (mre inactive) or log_sigma1/sigma2 (mre active).
  Va <- V[1 : k, 1 : k, drop=FALSE]
  Vsumm <- Matrix::summary(methods::as(Matrix::drop0(Va), "TsparseMatrix"))
  keep  <- Vsumm$i <= Vsumm$j
  a_list <- list(row = Vsumm$i[keep], col = Vsumm$j[keep], a = Vsumm$x[keep])

  b_mat  <- as.matrix(V[(k + 1) : (k + ptau), (k + 1) : (k + ptau), drop=FALSE])
  ab_mat <- as.matrix(V[1 : k, (k + 1) : (k + ptau), drop=FALSE])

  list(a = a_list, b = b_mat, ab = ab_mat,
      iname = iname,
      xname = c(xname, 'log(sigma)', if(use_mre) 'sigma2'))
}

## Claude Sonnet 5 2026-09-06
## mre added: optional fixed, user-specified "multiplier for random
## effects" template (same length as y). NULL (default) reproduces
## every existing formula and code path exactly -- p_wt is 1 (just
## log(sigma)), theta_full's layout is unchanged, and wtre is an
## all-ones vector wherever it's used. When supplied, TWO parameters
## are estimated JOINTLY with alpha/beta via the ordinary Newton step:
## sigma1 (log-linked, positive) and sigma2 (unconstrained), giving
## the random effect v_i~N(0,1) a weight of [sigma1*(1-mre)+sigma2*mre]
## at each observation -- see agqStep's header comment for the full
## derivation and for why this additive form, rather than one scale
## multiplicatively modifying another, is what makes both parameters
## well-identified. sigma2 starts at exactly 0 (not sigma.init) so the
## Newton step can discover its correct sign directly from the data --
## confirmed directly that starting it at a nonzero value biases the
## fit toward that sign, converging to a persistent wrong optimum when
## the truth has the opposite sign, regardless of sample size.
## The recommended usage for a Markov-1 model is mre=0 at each
## subject's first (anchor) observation and mre=1 thereafter, giving
## weight sigma1 at the anchor and sigma2 at every later observation --
## sigma2 of the same sign as sigma1 lets the random effect add to the
## lag-1-induced correlation, while a sigma2 of opposite sign lets it
## subtract from an over-induced one. mre need not vary in every
## cluster: a cluster with a single (anchor-only) observation
## contributes nothing about sigma1/sigma2 but causes no problem;
## they are identified as long as SOME clusters have within-cluster
## mre variation, checked directly below (mre constant at any single
## value, not just 0 or 1, is rejected).
ormrfit <- function(x, y, y2, k, intcens=0L, cluster, nc, initial, sigma.init=1.0,
                    offset=rep(0., n), wt=rep(1., n), penmat=matrix(0., p, p),
                    maxit=30L, maxit.outer=100L, maxit.mode=30L,
                    objtol=5e-4, gradtol=1e-3, paramtol=1e10,
                    tolsolve=.Machine$double.eps, minstepsize=1e-2, trace=FALSE,
                    link, iname, xname,
                    nAGQ=7L, nAGQ.grid=c(7L, 11L, 15L, 21L, 31L, 45L, 63L),
                    nAGQ.tol=1e-5, mre=NULL) {

  # Claude Sonnet 5 2026-08-30          entire function
  n <- length(y)
  p <- length(initial) - k
  deb <- Hmisc::Fdebug('orm.re.debug')
  deb(k); deb(length(initial)); deb(p); deb(dim(x)); deb(initial)

  storage.mode(x)       <- 'double'
  storage.mode(y)       <- 'integer'
  storage.mode(y2)      <- 'integer'
  storage.mode(k)       <- 'integer'
  storage.mode(p)       <- 'integer'
  storage.mode(cluster) <- 'integer'
  storage.mode(offset)  <- 'double'
  storage.mode(wt)      <- 'double'
  storage.mode(penmat)  <- 'double'
  storage.mode(link)    <- 'integer'

  ## Claude Sonnet 5 2026-09-06          10 lines
  use_mre <- length(mre) > 0
  if(use_mre) {
    if(length(mre) != n) stop('mre must have the same length as y')
    ## Claude Sonnet 5 2026-09-06
    ## sigma1/sigma2 are separately identifiable only if mre varies
    ## within the data at all -- constant mre (at 0, at 1, or at any
    ## other single value) leaves only one combination of sigma1/
    ## sigma2 identified, not both. diff(range(mre))==0 catches all
    ## three cases in one check (faster than length(unique(mre))<2).
    if(diff(range(mre)) == 0.0) {
      warning('mre has no variation (all values equal ', mre[1],
              '); sigma1 and sigma2 are not separately identifiable')
      return(list(fail=TRUE))
    }
    storage.mode(mre) <- 'double'
  }
  ## p_wt: total extra parameters beyond alpha/beta -- 2 (log_sigma1,
  ## sigma2) when mre is active, 1 (log(sigma)) otherwise. Always at
  ## least 1: there is always some random-effect scale being estimated.
  p_wt <- if(use_mre) 2L else 1L

  ## ia/ia2/sgn/ib/nb: pure function of y, y2, k -- computed once,
  ## shared by every ormll and ormeta call for the whole fit.
  idx <- .Fortran(F_ormidx, n, k, y, y2,
                  ia=integer(n), ia2=integer(n), sgn=numeric(n),
                  ib=integer(n), nb=integer(1), salloc=integer(1))
  if(idx$salloc != 0)
    stop('Censoring values encountered that are not handled (ormidx code ',
        idx$salloc, ')')
  nb <- idx$nb
  ## Claude Sonnet 5 2026-08-30          1 line
  ## Guard against R's 1:n gotcha (1:0 == c(1,0), not empty) if a
  ## dataset ever has nb=0 (no observations needing a second alpha).
  ib <- if(nb > 0) idx$ib[1 : nb] else integer(0)

  m <- function(v) max(abs(v))

  theta   <- initial
  alpha   <- theta[1 : k]
  beta    <- if(p > 0) theta[(k + 1) : (k + p)] else numeric(0)
  sigma   <- sigma.init
  ## Claude Sonnet 5 2026-09-06          6 lines
  ## sigma1 starts at sigma.init (needs some positive value to start
  ## the log-link from); sigma2 starts at exactly 0, NOT sigma.init.
  ## Confirmed directly: starting sigma2 at a nonzero value biased the
  ## fit toward whichever sign it started on, converging to a
  ## persistent wrong optimum when that sign didn't match the truth
  ## (robust to sample size, not fixed by more data). sigma2=0 has no
  ## such bias -- its own pseudo-covariate has no chain-rule scale
  ## factor (unlike sigma1's), so its gradient at exactly 0 is well-
  ## defined and lets the Newton step discover the correct sign
  ## directly from the data, verified on both positive- and negative-
  ## truth scenarios.
  if(use_mre) { sigma1 <- sigma.init; sigma2 <- 0 }
  gamma_c     <- rep(0e0, nc)
  xbeta   <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)
  base_lp <- offset + xbeta
  deb(beta); deb(xbeta); deb(base_lp)

  nAGQ.grid <- sort(unique(c(nAGQ, nAGQ.grid)))
  cur.nAGQ  <- nAGQ.grid[1]
  grid.pos  <- 1L
  ## Claude Sonnet 5 2026-09-04          1 line
  ## outer.iter resets to 0 at the top of every nAGQ escalation pass
  ## (needed so maxit.outer remains a per-nAGQ-level budget, matching
  ## the pre-SQUAREM loop's own semantics); total.oi accumulates
  ## across passes so the final reported $iter reflects the true total
  ## computational cost (raw onestep() calls) of the whole fit, not
  ## just the last nAGQ level's count.
  total.oi <- 0L

  repeat {   # nAGQ escalation: refit (warm-started) at each grid value
             # in turn until the log-likelihood stabilizes

    ## Claude Sonnet 5 2026-09-04          1 line
    ## Rebuilt fresh each pass through this outer repeat{} -- on the
    ## first pass from the function's initial theta/sigma, and on
    ## every subsequent pass (after an nAGQ escalation) from whatever
    ## theta/sigma the previous nAGQ level converged to, since those
    ## are updated as plain variables right after the accelerated
    ## loop below breaks.
    ## Claude Sonnet 5 2026-09-06          6 lines
    ## theta_full's layout is c(alpha, beta, [log_sigma1, sigma2]) when
    ## mre is active, or c(alpha, beta, log(sigma)) otherwise -- these
    ## are mutually exclusive, not additive: when mre is active there
    ## is no separate "sigma" at all (v_i's prior is fixed at 1), so
    ## sigma1/sigma2 together replace it rather than sitting alongside
    ## it. log_sigma1 keeps sigma1 positive; sigma2 is stored directly
    ## (unconstrained, no transform).
    theta_full <- if(use_mre) c(theta, log(sigma1), sigma2) else c(theta, log(sigma))

    ## Claude Sonnet 5 2026-09-04          entire onestep() function
    ## One full outer "fixed-point" iteration -- per-cluster
    ## mode-finding, the AGQ-reweighted Newton step for (alpha, beta,
    ## [sigma1, sigma2]) with its own step-halving, and (when mre is
    ## inactive) the sigma update -- factored out unchanged from the
    ## previous plain loop so it can be called repeatedly per SQUAREM
    ## cycle below. grad/delta in its return value are evaluated AT
    ## theta_full going IN (matching the pre-existing semantics
    ## exactly: gradient before the step, size of that step), not at
    ## the result.
    ## Claude Sonnet 5 2026-09-06          revision note
    ## sigma (as tau=log(sigma), the mre-INACTIVE case) is part of the
    ## SAME joint Newton step as alpha/beta, rather than a separate 1-D
    ## profile search. Motivation, confirmed directly during an earlier
    ## (now superseded) design where the random effect was weighted by
    ## a single re_mult coefficient: a separate-search scheme let a
    ## step improve the alpha/beta/re_mult surrogate while the
    ## SUBSEQUENT, independently-chosen sigma update actually made the
    ## true AGQ log-likelihood WORSE overall -- verified by direct
    ## comparison (the true parameter values gave a clearly higher
    ## log-likelihood than the point the old scheme was converging
    ## toward) in a scenario with a true re_mult>1 (a negative net
    ## weight). Folding tau into the same step-halving check closes
    ## that gap: any joint step is now only accepted if it doesn't
    ## degrade the SAME surrogate tau is now part of. This same
    ## motivation is why sigma1/sigma2 (the mre-ACTIVE case) were
    ## designed from the start to be part of this joint step too,
    ## rather than needing a similar retrofit.
    ##
    ## tau's own gradient/curvature (accum_tau_grad/accum_tau_curv in
    ## agqStep) are computed directly from the AGQ node weights and
    ## values, not via ormll -- tau affects only the random-effect
    ## PRIOR term, not eta, so it isn't a pseudo-covariate the way
    ## sigma1/sigma2 are. Verified against finite differences of the
    ## EM-surrogate objective (exact agreement). Under this fixed-mode
    ## surrogate, tau's cross-terms with alpha/beta are exactly zero
    ## (the data log-lik term doesn't depend on sigma and the prior
    ## term doesn't depend on alpha/beta, once gamma_node_c is held
    ## fixed) -- so the extended Hessian is block-diagonal. This does
    ## not capture the TRUE curvature coupling between sigma and the
    ## other parameters (which flows through the mode gamma_c/kappa_c
    ## itself, not through this fixed-
    ## mode surrogate) -- but it does put sigma under the same trial-
    ## and-reject step-halving discipline as everything else, which is
    ## what the failure above actually needed.
    onestep <- function(theta_full, gamma_c, oi, role) {
      alpha   <- theta_full[1 : k]
      beta    <- if(p > 0) theta_full[(k + 1) : (k + p)] else numeric(0)
      ## Claude Sonnet 5 2026-09-06          8 lines
      ## Mutually exclusive extraction: when mre is active, sigma1/
      ## sigma2 are the last two positions (no separate tau at all,
      ## and the mode-finding/AGQ prior SD for v_i is fixed at 1, not
      ## estimated); otherwise the last position is tau=log(sigma) as
      ## before.
      if(use_mre) {
        sigma1      <- exp(theta_full[k + p + 1])
        sigma2      <- theta_full[k + p + 2]
        sigma1_prev <- sigma1   ## Claude Sonnet 5 2026-09-06: kept separately from `sigma`
        sigma       <- 1        ## below, since `sigma` is overloaded as v_i's FIXED prior SD
      } else {                  ## here, not the previous iteration's estimated scale.
        sigma1 <- NULL; sigma2 <- NULL
        sigma  <- exp(theta_full[k + p + p_wt])
      }
      xbeta   <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)
      base_lp <- offset + xbeta
      wtre    <- if(use_mre) sigma1 * (1 - mre) + sigma2 * mre else rep(1e0, n)

      mf <- clusterModeFind(alpha, beta, sigma, gamma_c, base_lp,
                            wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                            cluster, nc, link, k, n, wtre=wtre, maxit=maxit.mode)
      if(mf$fail) {
        if(trace) message('clusterModeFind failed (ormeta salloc=', mf$code,
                          ') at nAGQ=', cur.nAGQ, ' outer.iter=', oi,
                          ' sigma=', sigma)
        return(list(fail=TRUE))
      }
      gamma_c <- mf$gamma_c

      ## Claude Sonnet 5 2026-08-30          3 lines
      ## Print the mode-finding/curvature state as soon as it's known,
      ## not just after a successful Newton step -- a failure further
      ## down (agqStep, the Hessian solve) would otherwise leave a
      ## trace=2 run with no output at all for that outer iteration.
      if(trace > 1)
        cat('  [nAGQ:', cur.nAGQ, ' outer.iter:', oi, ' sigma:', sigma,
            ' kappa_c range:', range(mf$kappa_c),
            ' step_c range:', range(sqrt(2 / mf$kappa_c)), ']\n')
      if(use_mre && trace > 1) cat('  [sigma1:', sigma1, ' sigma2:', sigma2, ']\n')

      ag <- agqStep(alpha, beta, sigma, gamma_c, mf$kappa_c, cur.nAGQ,
                   base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                   cluster, nc, link, k, p, n, penmat, penhess=1L,
                   intcens=intcens, score=TRUE, mre=mre, sigma1=sigma1, sigma2=sigma2)
      if(ag$fail) {
        if(trace) message('agqStep failed (code=', ag$code,
                          ') at nAGQ=', cur.nAGQ, ' outer.iter=', oi,
                          ' sigma=', sigma)
        return(list(fail=TRUE))
      }
      if(trace > 1)
        cat('  [loglik:', ag$loglik, ' grad:', ag$grad, ']\n')

      ## Claude Sonnet 5 2026-08-30          2 lines
      ## infoMxop already natively supports a sparse triplet list for
      ## `a` (its own documented format for ormfit's intcens=1 case) --
      ## just construct the right shape depending on intcens.
      a_for_hess <- if(intcens == 1L) list(row=ag$row, col=ag$col, a=ag$ai) else ag$ha
      hess <- infoMxop(list(a=a_for_hess, b=ag$hb, ab=ag$hab))
      ## Claude Sonnet 5 2026-09-07          revision note
      ## Reverted: sigma (tau) is NO LONGER folded into this joint step
      ## when mre is inactive. That fold was added to fix a genuine
      ## problem in the mre-active (sigma1/sigma2) case, but confirmed
      ## directly to REGRESS the far more common plain (mre-inactive)
      ## case on real data: a single linearized Newton step for
      ## log(sigma) can systematically undershoot when sigma's own
      ## curvature is asymmetric, producing a slow, monotonic,
      ## never-converging drift (sigma steadily decreasing for 100
      ## iterations without stabilizing, where the original separate-
      ## profile-search design converged in 3). hess/ag$grad are used
      ## directly here -- sigma1/sigma2 are ALREADY part of them (via
      ## p_star=p+2 in agqStep) when mre is active, and sigma is
      ## handled by its own separate 1-D search (sigmaObjective) below
      ## when it is not.
      full_grad <- ag$grad

      ## Claude Sonnet 5 2026-08-30          10 lines
      ## Levenberg-Marquardt-style ridge fallback, matching ormfit's
      ## own established formula exactly (hess + lambda*diag(hess), NOT
      ## a plain identity -- hess is on the natural +LL scale, negative
      ## definite at a good point, so scaling by its own diagonal is
      ## what correctly moves it further from singular as lambda grows;
      ## a plain +lambda*I would push the wrong direction on this
      ## sign convention). A near-singular accumulated Hessian is
      ## plausible in harder regimes (large sigma, sharp per-cluster
      ## curvature concentrating almost all AGQ weight on one node) --
      ## retry with increasing damping before concluding real failure.
      ## Claude Sonnet 5 2026-09-06          3 lines
      ## Ridge scale floored at 1: hess + lambda*diag(hess) does
      ## nothing at any position where hess's own diagonal is exactly
      ## zero (lambda*0=0 regardless of lambda) -- a real, reachable
      ## failure mode confirmed directly during the earlier (wl,rho)
      ## work, retained here as a general safeguard.
      lambda <- 0
      repeat {
        hess_try <- if(lambda == 0) hess else
                      hess + lambda * Matrix::Diagonal(x=pmax(abs(Matrix::diag(hess)), 1))
        delta <- try(Matrix::solve(hess_try, full_grad, tol=tolsolve), silent=TRUE)
        if(! inherits(delta, 'try-error')) break
        lambda <- if(lambda == 0) 1e-6 else lambda * 10
        if(lambda > 1e6) {
          if(trace) {
            message('singular Hessian matrix in ormrfit even after ridge damping ',
                    'at nAGQ=', cur.nAGQ, ' outer.iter=', oi, ' sigma=', sigma)
            ev <- try(eigen(as.matrix(hess), only.values=TRUE)$values, silent=TRUE)
            if(! inherits(ev, 'try-error'))
              cat('  [hess eigenvalues:', ev, ']\n')
            cat('  [hess diag:', Matrix::diag(hess), ']\n')
            ## Claude Sonnet 5 2026-08-30          6 lines
            ## Distinguish "beta's whole row/column is zero" (points to
            ## x not reaching the accumulated ormll calls) from "just
            ## the diagonal happens to vanish" (could be a genuine
            ## envelope-theorem-style cancellation at this iteration).
            cat('  [full hess:\n'); print(as.matrix(hess)); cat('  ]\n')
            cat('  [beta:', beta, ' range(x):', range(x), ' table(x):',
                table(x), ']\n')
          }
          return(list(fail=TRUE))
        }
      }
      if(trace > 1 && lambda > 0)
        cat('  [ridge damping used, lambda=', lambda, ']\n')

      ## Step-halving surrogate: an EM/MM-style "M-step" objective --
      ## the node posterior weights (ag$rho) and the mode/curvature
      ## (gamma_c, ag$kappa_c / ag$step_c / ag$gh) are held FIXED at their
      ## values from the current outer iteration while trialing a new
      ## theta, giving sum_m rho_m * clusterLogLik(theta, node_m),
      ## which shares a gradient with the true AGQ log-likelihood at
      ## the current theta and is cheap (ormeta-only) to evaluate
      ## repeatedly. This is NOT the true log-sum-exp AGQ likelihood
      ## away from the current theta -- that gets recomputed properly
      ## at the *start* of every outer iteration via agqStep, so any
      ## gap between the surrogate and the true objective cannot
      ## accumulate silently across outer iterations.
      ## Claude Sonnet 5 2026-09-07
      ## Trials (alpha,beta) alone when mre is inactive (sigma is NOT
      ## part of this step, see revision note above -- no prior-density
      ## term needed here, matching the original pre-fold design), or
      ## (alpha,beta,log_sigma1,sigma2) jointly when mre is active,
      ## INCLUDING the prior-density term (sigma1/sigma2 genuinely are
      ## part of this step in that case, exactly as validated for the
      ## negative-weight scenario).
      theta_star <- if(use_mre) theta_full else theta_full[1 : (k + p)]
      surrogate <- function(theta_t) {
        alpha_t <- theta_t[1 : k]
        if(k > 1 && any(diff(alpha_t) >= 0)) return(Inf)
        beta_t    <- if(p > 0) theta_t[(k + 1) : (k + p)] else numeric(0)
        xbeta_t   <- if(p > 0) as.vector(x %*% beta_t) else rep(0e0, n)
        lp_base_t <- offset + xbeta_t
        if(use_mre) {
          sigma1_t <- exp(theta_t[k + p + 1])
          sigma2_t <- theta_t[k + p + 2]
          wtre_t   <- sigma1_t * (1 - mre) + sigma2_t * mre
        } else wtre_t <- rep(1e0, n)
        tot <- 0
        for(mm in 1 : cur.nAGQ) {
          gamma_node_c <- gamma_c + ag$step_c * ag$gh$nodes[mm]
          gamma_node   <- wtre_t * gamma_node_c[cluster]
          lp_t     <- lp_base_t + gamma_node
          ## Claude Sonnet 5 2026-09-08          2 lines
          ## Same non-finite-argument guard as clusterModeFind/agqStep,
          ## returning Inf here to match this function's own existing
          ## graceful-rejection convention (salloc!=0 -> Inf), so the
          ## surrogate step-halving loop above simply treats this trial
          ## step size as rejected and halves further, rather than
          ## crashing on an uncatchable .Fortran() error.
          if(any(! is.finite(lp_t))) return(Inf)
          wtt <- .Fortran(F_ormeta, n, k, link, alpha_t, lp1=lp_t, lp2=lp_t, wt,
                          idx$ia, idx$ia2, idx$sgn, ib, nb,
                          logd=numeric(n), g=numeric(n), h=numeric(n),
                          s1=numeric(n), s2=numeric(n),
                          debug=0L, salloc=integer(1))
          if(wtt$salloc != 0) return(Inf)
          tot <- tot + sum(ag$rho[cluster, mm] * wtt$logd)
          if(use_mre)
            tot <- tot + sum(ag$rho[, mm] * dnorm(gamma_node_c, mean=0, sd=1, log=TRUE))
        }
        -2 * tot
      }

      objf      <- surrogate(theta_star)
      step_size <- 1.0
      repeat {
        theta_star_new <- theta_star - step_size * delta
        objfnew   <- surrogate(theta_star_new)
        if(trace > 1)
          cat('    [step_size=', step_size, ' objf=', objf, ' objfnew=', objfnew, ']\n')
        if(! is.finite(objfnew) || objfnew > objf + objtol / 10) {
          step_size <- step_size / 2
          if(step_size < minstepsize) {
            message('step size reduced below minstepsize in ormrfit ',
                    'without improving the objective')
            return(list(fail=TRUE))
          }
        } else break
      }

      alpha_new <- theta_star_new[1 : k]
      beta_new  <- if(p > 0) theta_star_new[(k + 1) : (k + p)] else numeric(0)
      xbeta_new <- if(p > 0) as.vector(x %*% beta_new) else rep(0e0, n)
      base_lp_new <- offset + xbeta_new

      ## Claude Sonnet 5 2026-09-07          entire block
      ## Re-added: when mre is inactive, sigma gets its own separate
      ## 1-D profile search (exactly as before the ill-fated joint-step
      ## fold) -- alpha/beta's step above never touched sigma, and
      ## sigma1/sigma2 (when mre IS active) were already updated as
      ## part of theta_star_new above, needing no further search here.
      if(use_mre) {
        sigma1_new <- exp(theta_star_new[k + p + 1])
        sigma2_new <- theta_star_new[k + p + 2]
        theta_full_new <- theta_star_new
      } else {
        br   <- log(sigma) + c(-2.5, 2.5)
        sopt <- optimize(sigmaObjective, interval=br, maximum=TRUE,
                         alpha=alpha_new, beta=beta_new, gamma_c_init=gamma_c, base_lp=base_lp_new,
                         wt=wt, ia=idx$ia, ia2=idx$ia2, sgn=idx$sgn, ib=ib, nb=nb,
                         cluster=cluster, nc=nc, link=link, k=k, n=n,
                         nAGQ=cur.nAGQ, maxit.mode=maxit.mode, tol=1e-4)
        sigma_new_local <- exp(sopt$maximum)
        theta_full_new <- c(theta_star_new, log(sigma_new_local))
      }

      ## Claude Sonnet 5 2026-09-06          4 lines
      ## sigma_new tracks sigma1 (the constrained, log-linked scale)
      ## when mre is active, matching the ORIGINAL boundary-oscillation
      ## concern this check exists for -- sigma2 is unconstrained and
      ## has no analogous near-zero-boundary degeneracy to guard against.
      sigma_new <- if(use_mre) exp(theta_full_new[k + p + 1]) else exp(theta_full_new[k + p + p_wt])
      sigma_prev_for_check <- if(use_mre) sigma1_prev else sigma
      sigma.change <- abs(sigma_new - sigma_prev_for_check) / max(sigma_prev_for_check, 1e-6)
      ## Claude Sonnet 5 2026-08-30          6 lines
      ## Near a sigma-to-zero boundary, sigma can oscillate in pure
      ## numerical noise (e.g. between 1e-7 and 1e-6) while the
      ## RELATIVE change above stays permanently large, since
      ## max(sigma,1e-6) floors the denominator right where sigma
      ## itself is already negligible -- this kept the outer loop
      ## running to maxit.outer even when grad/delta were already at
      ## 1e-13/1e-14 (genuinely converged). Also accept convergence
      ## when the ABSOLUTE change is below that same floor, which only
      ## ever matters once sigma is already numerically at the
      ## boundary -- harmless for any normal, well-identified sigma.
      sigma.converged <- sigma.change < objtol || abs(sigma_new - sigma) < 1e-6

      ## Claude Sonnet 5 2026-09-04          5 lines
      ## Moved out of an unconditional print: onestep() is now called
      ## 2-3 times per SQUAREM cycle (s1, s2, and an optional
      ## stabilization trial that may be REJECTED), and printing every
      ## call unlabeled made a rejected trial's wild intermediate
      ## values (confirmed: sigma swinging to 0.53 with delta=33.7 on
      ## a real trace) look like an accepted, alarming divergence. The
      ## real per-cycle summary (accepted values only) is printed by
      ## the cycle loop below; this detailed, explicitly-labeled
      ## per-call trace is opt-in via trace>1.
      if(trace > 1)
        cat('  [', role, ' nAGQ:', cur.nAGQ, ' outer.iter:', oi,
            ' sigma:', format(sigma_new, nsmall=4),
            ' max|grad|:', m(full_grad), ' max|delta|:', m(delta), ']\n')

      list(fail=FALSE, theta_full=theta_full_new,
          gamma_c=gamma_c, grad=full_grad, delta=delta,
          sigma.converged=sigma.converged)
    }



    ## Claude Sonnet 5 2026-09-04          entire loop replacement
    ## SQUAREM acceleration (Varadhan & Roland 2008, "SqS3" scheme)
    ## for the outer loop. Motivation: on a real n=20 test case this
    ## loop showed clean, purely geometric convergence (max|delta|
    ## shrinking by a near-constant ~0.86 per iteration for dozens of
    ## iterations) -- exactly the situation this kind of extrapolation
    ## is designed for, and needed ~57 plain iterations to reach
    ## paramtol where a handful of SQUAREM cycles should suffice.
    ## Mechanism: two plain onestep() calls give theta0, theta1,
    ## theta2; r=theta1-theta0 and v=(theta2-theta1)-r estimate the
    ## local geometric contraction, giving a closed-form extrapolated
    ## point theta_sq. That point is never accepted on faith -- it is
    ## itself run through one more onestep() ("stabilization"), and
    ## only kept if that step's own delta is at least as small as the
    ## plain iterate's (s2's). Any invalid extrapolation (non-finite,
    ## non-descending alpha, degenerate v, or the stabilization step
    ## itself failing) falls back to the plain s2 -- so this can only
    ## ever match or beat plain iteration, never do meaningfully worse.
    ## outer.iter counts total onestep() calls (matching the previous
    ## per-plain-iteration meaning of maxit.outer as a cost budget) and
    ## is capped from being exceeded by a stabilization trial.
    outer.iter <- 0L
    hit.maxit  <- FALSE
    stepmax    <- 1
    repeat {
      s1 <- onestep(theta_full, gamma_c, outer.iter + 1L, 's1')
      outer.iter <- outer.iter + 1L
      if(s1$fail) return(list(fail=TRUE))

      if(outer.iter >= maxit.outer) {
        accepted <- s1
        hit.maxit <- TRUE
      } else {
        s2 <- onestep(s1$theta_full, s1$gamma_c, outer.iter + 1L, 's2')
        outer.iter <- outer.iter + 1L
        if(s2$fail) return(list(fail=TRUE))

        accepted <- s2   # default: plain second iterate
        r  <- s1$theta_full - theta_full
        v  <- (s2$theta_full - s1$theta_full) - r
        nr <- sqrt(sum(r * r)); nv <- sqrt(sum(v * v))
        if(nv > 1e-12 * max(1, nr) && outer.iter < maxit.outer) {
          ## Claude Sonnet 5 2026-09-04          12 lines
          ## Standard SQUAREM safeguard (Varadhan & Roland 2008):
          ## bound the extrapolation's step length rather than only
          ## detecting a bad result after paying for a full onestep()
          ## trial. |step_len|=nr/nv is unbounded when the local
          ## trajectory isn't yet in a clean geometric regime (e.g.
          ## early cycles, or when different components of theta_full
          ## converge at different rates) -- this is what produced the
          ## earlier wild trial (sigma swinging to 0.53, delta=33.7).
          ## stepmax starts conservative and grows 4x after each
          ## ACCEPTED extrapolation (confidence the local geometric
          ## model is trustworthy), shrinking back 4x after a
          ## rejection -- so wasted trials become rarer as the fit
          ## settles into its true linear-convergence regime.
          step_len <- max(-nr / nv, -stepmax)
          theta_sq <- theta_full - 2 * step_len * r + step_len^2 * v
          alpha_sq <- theta_sq[1 : k]
          valid <- all(is.finite(theta_sq)) && (k == 1 || all(diff(alpha_sq) < 0))
          if(valid) {
            s3 <- onestep(theta_sq, s2$gamma_c, outer.iter + 1L, 'squarem-trial')
            outer.iter <- outer.iter + 1L
            squarem.accepted <- ! s3$fail && m(s3$delta) <= m(s2$delta)
            if(squarem.accepted) accepted <- s3
            stepmax <- if(squarem.accepted) stepmax * 4 else max(1, stepmax / 4)
            if(trace > 1)
              cat('  [squarem-trial', if(squarem.accepted) 'ACCEPTED' else 'rejected',
                  '-- plain m(delta)=', m(s2$delta),
                  ' trial m(delta)=', if(s3$fail) NA else m(s3$delta),
                  ' step_len=', step_len, ' stepmax now', stepmax, ']\n')
          }
        }
        if(outer.iter >= maxit.outer) hit.maxit <- TRUE
      }

      theta_full <- accepted$theta_full
      gamma_c    <- accepted$gamma_c

      ## Claude Sonnet 5 2026-09-04          3 lines
      ## The real, "official" per-cycle summary -- accepted values
      ## only. onestep() itself no longer prints unconditionally (see
      ## its own comment above): a rejected trial's wild intermediate
      ## values were previously indistinguishable from genuine
      ## progress in a trace=1 run.
      if(trace) {
        sigma_trace <- if(use_mre) exp(theta_full[k + p + 1]) else exp(theta_full[k + p + p_wt])
        cat('nAGQ:', cur.nAGQ, ' outer iter:', outer.iter,
            ' sigma:', format(sigma_trace, nsmall=4),
            if(use_mre) c(' sigma2:', format(theta_full[k + p + 2], nsmall=4)),
            ' max|grad|:', m(accepted$grad), ' max|delta|:', m(accepted$delta), '\n')
      }

      ## Claude Sonnet 5 2026-08-30          4 lines
      ## Unlike orm.fit's Newton-Raphson (quadratic convergence, so a
      ## loose n-scaled gradtol never actually binds -- it overshoots
      ## to near machine precision regardless), this outer loop is a
      ## block-alternating scheme (mode-find / AGQ-reweighted Newton /
      ## sigma profile) with only linear convergence near the optimum.
      ## An n-scaled gradient tolerance here is nearly a no-op at
      ## realistic sample sizes and lets the loop stop well short of
      ## genuine convergence, relying on sigma.change alone (which
      ## doesn't verify alpha/beta have also stabilized) to signal
      ## "done". Use gradtol as an absolute tolerance instead.
      if(m(accepted$grad) < gradtol && m(accepted$delta) < paramtol &&
        accepted$sigma.converged) break

      ## Claude Sonnet 5 2026-09-03          6 lines
      ## Matches ormfit's own established convention exactly (see its
      ## "Reached maxit iterations without convergence" checks in both
      ## its NR and L-M branches): exhausting the outer loop without the
      ## convergence break ever firing is treated as a hard failure, not
      ## silently returned as if converged. Previously this fell through
      ## into the nAGQ escalation check and the final return with
      ## fail=FALSE regardless, giving no signal that the fit was less
      ## precise than it looked.
      if(hit.maxit) {
        sigma_msg <- if(use_mre) exp(theta_full[k + p + 1]) else exp(theta_full[k + p + p_wt])
        msg <- paste('Reached', maxit.outer, 'outer iterations without convergence at nAGQ=',
                    cur.nAGQ, '\nsigma:', sigma_msg,
                    ' Max |gradient|:', m(accepted$grad),
                    ' Max |change in parameters|:', m(accepted$delta))
        message(msg)
        return(list(fail=TRUE))
      }
    }

    total.oi <- total.oi + outer.iter

    alpha   <- theta_full[1 : k]
    beta    <- if(p > 0) theta_full[(k + 1) : (k + p)] else numeric(0)
    if(use_mre) {
      sigma1 <- exp(theta_full[k + p + 1])
      sigma2 <- theta_full[k + p + 2]
      sigma  <- 1
    } else {
      sigma  <- exp(theta_full[k + p + p_wt])
    }
    theta   <- theta_full[1 : (k + p)]
    xbeta   <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)
    base_lp <- offset + xbeta

    ## Escalation check: does the next larger nAGQ change the
    ## log-likelihood by more than nAGQ.tol? (cheap, ormeta-only,
    ## reusing the mode/theta just converged to as a warm start)
    wtre.cur <- if(use_mre) sigma1 * (1 - mre) + sigma2 * mre else rep(1e0, n)
    mf.cur <- clusterModeFind(alpha, beta, sigma, gamma_c, base_lp,
                              wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                              cluster, nc, link, k, n, wtre=wtre.cur, maxit=maxit.mode)
    ## Claude Sonnet 5 2026-09-08          entire block
    ## This escalation check is a cheap, OPTIONAL verification step run
    ## AFTER the outer loop has already converged at cur.nAGQ -- if
    ## clusterModeFind or either agqStep call below fails (e.g. a
    ## wider-spread AGQ node under the larger candidate nAGQ hits a
    ## genuinely invalid category probability -- confirmed to occur
    ## under cloglog), there is no meaningful rel.change to compute at
    ## all. agqStep's own failure return, list(fail=TRUE, code=...),
    ## has no $loglik element, so $loglik on it is NULL -- and
    ## abs(NULL - ll.cur) silently collapses to numeric(0), which then
    ## crashes the if() below with exactly "argument is of length
    ## zero" rather than any error naming the real cause. Rather than
    ## let that happen, explicitly check each call's own $fail and
    ## treat an inability to perform this check as "cannot verify
    ## stability at a larger nAGQ" -- simply keep the already-converged
    ## result at cur.nAGQ, exactly as if the escalation grid had
    ## reached its ceiling.
    if(mf.cur$fail) {
      if(trace) message('nAGQ escalation check skipped (clusterModeFind failed at cur.nAGQ=',
                        cur.nAGQ, ') -- keeping the already-converged result')
      break
    }
    agq.cur <- agqStep(alpha, beta, sigma, mf.cur$gamma_c, mf.cur$kappa_c, cur.nAGQ,
                       base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                       cluster, nc, link, k, p, n, penmat, score=FALSE,
                       mre=mre, sigma1=sigma1, sigma2=sigma2)
    if(agq.cur$fail) {
      if(trace) message('nAGQ escalation check skipped (agqStep failed at cur.nAGQ=',
                        cur.nAGQ, ') -- keeping the already-converged result')
      break
    }
    ll.cur <- agq.cur$loglik

    if(grid.pos == length(nAGQ.grid)) {
      if(trace) message('nAGQ reached grid ceiling (', cur.nAGQ,
                        ') without full convergence -- results may still be changing')
      break
    }
    next.nAGQ <- nAGQ.grid[grid.pos + 1]
    agq.next <- agqStep(alpha, beta, sigma, mf.cur$gamma_c, mf.cur$kappa_c, next.nAGQ,
                        base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                        cluster, nc, link, k, p, n, penmat, score=FALSE,
                        mre=mre, sigma1=sigma1, sigma2=sigma2)
    if(agq.next$fail) {
      if(trace) message('nAGQ escalation check inconclusive (agqStep failed at next.nAGQ=',
                        next.nAGQ, ') -- keeping the already-converged result at cur.nAGQ=',
                        cur.nAGQ)
      break
    }
    ll.next   <- agq.next$loglik
    rel.change <- abs(ll.next - ll.cur) / max(abs(ll.cur), 1e-6)
    ## Claude Sonnet 5 2026-08-30          5 lines
    ## Bug fix: previously set cur.nAGQ <- next.nAGQ here before
    ## breaking -- but the whole point of this check passing is that
    ## the CURRENT (smaller) nAGQ was already adequate; theta/sigma
    ## were only ever driven to convergence under that value, not
    ## under next.nAGQ, so reporting/using next.nAGQ here left the
    ## final gradient reflecting a mismatch between the converged
    ## theta and a different (if very similar) objective it was never
    ## actually optimized against. Stay at the current, converged nAGQ.
    if(rel.change < nAGQ.tol) break

    grid.pos <- grid.pos + 1L
    cur.nAGQ <- next.nAGQ
    gamma_c      <- mf.cur$gamma_c   # warm start into the next nAGQ's refit
  }

  ## Final accumulation at converged (alpha, beta, sigma, nAGQ) for
  ## the returned information matrix. No special "include_wt_cols"
  ## flag needed here (unlike the earlier wl/rho design): sigma1/
  ## sigma2's Hessian contribution is already correctly included in
  ## ag.final's own hb/hab, since it was jointly Newton-updated
  ## throughout.
  wtre.final <- if(use_mre) sigma1 * (1 - mre) + sigma2 * mre else rep(1e0, n)
  mf.final <- clusterModeFind(alpha, beta, sigma, gamma_c, base_lp,
                              wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                              cluster, nc, link, k, n, wtre=wtre.final, maxit=maxit.mode)
  ag.final <- agqStep(alpha, beta, sigma, mf.final$gamma_c, mf.final$kappa_c, cur.nAGQ,
                      base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                      cluster, nc, link, k, p, n, penmat,
                      intcens=intcens, score=TRUE, mre=mre, sigma1=sigma1, sigma2=sigma2)

  info <- sparseInfoMatrix(alpha, beta, sigma, mf.final$gamma_c, mf.final$kappa_c,
                           base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                           cluster, nc, link, k, p, n,
                           ag.final$ha, ag.final$hb, ag.final$hab,
                           ag.final$row, ag.final$col, ag.final$ai, ag.final$ne,
                           intcens, cur.nAGQ, iname, xname,
                           mre=mre, sigma1=sigma1, sigma2=sigma2)

  ## Claude Sonnet 5 2026-08-30
  ## Field names/shape aligned with ormfit's own return object as
  ## closely as possible, so orm.fit()'s post-processing (stats block,
  ## retlist construction) can treat both workers' output uniformly:
  ## coef/loglik/u/info/iter/fail match ormfit's naming exactly
  ## (grad -> u). ncluster/nAGQ/sigma/gamma have no ormfit analog and
  ## are orm.fit()'s cue to attach them only when cluster is present.
  ##
  ## Two fields ormfit provides that this does NOT yet meaningfully
  ## compute, flagged explicitly rather than silently omitted:
  ## - score: ormfit's score-test statistic is only defined at the
  ##   FIRST Newton iteration with beta==0 -- a precondition that
  ##   never holds here, since ormrfit is always warm-started from an
  ##   already-fit (generally nonzero-beta) unclustered model. Kept as
  ##   NA, which orm.fit's existing stats block already treats as a
  ##   normal, handled case (same as its own initial.there branch).
  ## - lpe: ormfit's per-observation "probability of the observed/
  ##   censored outcome" (used only for the anycens ESS adjustment in
  ##   orm.fit's stats block). No per-observation analog is computed
  ##   here (ormrfit tracks per-CLUSTER, AGQ-marginalized quantities,
  ##   not a per-observation marginal probability, which isn't as
  ##   clean a concept once observations share a random effect). Kept
  ##   as NA so anycens-adjusted ESS silently comes out NA for
  ##   clustered censored fits rather than using a wrong value -- a
  ##   known, deliberately-scoped-out gap, not an oversight.
  ##
  ## Claude Sonnet 5 2026-09-06          3 lines
  ## sigma1/sigma2 reported directly only when mre was used; NULL
  ## otherwise (orm.fit's own dispatch treats their presence as the
  ## cue to attach them to the returned fit object, exactly like sigma/
  ## gamma already are). `sigma` itself is NOT meaningful when mre is
  ## active (it was fixed at 1 throughout, v_i's own prior SD) -- NULL
  ## in that case rather than reporting the uninformative fixed value.
  list(coef    = theta,
      sigma    = if(use_mre) NULL else sigma,
      sigma1   = if(use_mre) sigma1 else NULL,
      sigma2   = if(use_mre) sigma2 else NULL,
      gamma    = mf.final$gamma_c,
      nAGQ     = cur.nAGQ,
      ncluster = nc,
      loglik   = -2 * ag.final$loglik,
      u        = ag.final$grad,
      dmax     = max(abs(ag.final$grad)),
      score    = NA_real_,
      lpe      = rep(NA_real_, n),
      mscore   = NULL,
      info     = info,
      iter     = total.oi,
      fail     = FALSE)
}
