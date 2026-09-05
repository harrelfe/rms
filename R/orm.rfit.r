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
## Claude Sonnet 5 2026-08-30          entire function
clusterModeFind <- function(alpha, beta, sigma, gamma_c_init, base_lp,
                            wt, ia, ia2, sgn, ib, nb, cluster, nc,
                            link, k, n, maxit=30L, tol=1e-8) {
  gamma_c <- gamma_c_init
  for(it in 1 : maxit) {
    gamma_obs <- gamma_c[cluster]
    lp  <- base_lp + gamma_obs
    w   <- .Fortran(F_ormeta, n, k, link, alpha, lp1=lp, lp2=lp, wt,
                    ia, ia2, sgn, ib, nb,
                    logd=numeric(n), g=numeric(n), h=numeric(n),
                    s1=numeric(n), s2=numeric(n),
                    debug=0L, salloc=integer(1))
    if(w$salloc != 0) return(list(fail=TRUE, code=w$salloc))
    Sg    <- as.vector(rowsum(w$g, cluster, reorder=TRUE))
    Sh    <- as.vector(rowsum(w$h, cluster, reorder=TRUE))
    grad  <- Sg - gamma_c / sigma^2
    curv  <- Sh - 1 / sigma^2        # negative at a proper maximum
    curv  <- pmin(curv, -1e-8)       # guard against numerical breakdown
    step  <- grad / curv
    gamma_c_new <- gamma_c - step
    if(max(abs(gamma_c_new - gamma_c)) < tol) { gamma_c <- gamma_c_new; break }
    gamma_c <- gamma_c_new
  }
  list(fail=FALSE, gamma_c=gamma_c, kappa_c=-curv, iter=it)
}

## AGQ log-likelihood (and, if score=TRUE, the accumulated score/Hessian
## for alpha & beta) at given (alpha, beta, sigma), using a pre-computed
## mode (gamma_c, kappa_c). ormeta supplies the cheap per-node cluster
## log-likelihoods; ormll (called once per node with reweighted wt and
## node-shifted offset) supplies the score/Hessian, reusing its
## existing, verified accumulation machinery unchanged. The penalty
## matrix is passed as zero to each per-node ormll call and applied
## exactly once afterward, to avoid counting it nAGQ times.
##
## base_lp here means offset + x %*% beta (full, as in clusterModeFind)
## when building ormeta's lp1/lp2; when building ormll's own `offset`
## argument, x %*% beta must be subtracted back out first, since ormll
## adds it itself internally -- see node_off below.
## Claude Sonnet 5 2026-08-30          entire function
agqStep <- function(alpha, beta, sigma, gamma_c, kappa_c, nAGQ,
                    base_lp, x, wt, ia, ia2, sgn, ib, nb,
                    cluster, nc, link, k, p, n, penmat, penhess=1L,
                    intcens=0L, score=TRUE) {
  gh     <- gauss_hermite_quad(nAGQ)
  step_c <- sqrt(2 / kappa_c)
  zeromat <- matrix(0e0, p, p)

  lt <- matrix(0e0, nc, nAGQ)
  for(m in 1 : nAGQ) {
    gamma_node_c <- gamma_c + step_c * gh$nodes[m]
    gamma_node   <- gamma_node_c[cluster]
    lp       <- base_lp + gamma_node
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
  rho       <- exp(lt - mrow)
  denom     <- rowSums(rho)
  loglik_c  <- log(step_c) + mrow + log(denom)   # natural-log scale
  rho       <- rho / denom                        # nc x nAGQ, rows sum to 1
  total_ll  <- sum(loglik_c)

  if(! score)
    return(list(fail=FALSE, loglik=total_ll, rho=rho, gh=gh, step_c=step_c,
               gamma_c=gamma_c, kappa_c=kappa_c))

  xbeta <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)

  nai <- if(intcens == 1L) 1000000L else 1L

  accum_grad <- numeric(k + p)
  accum_ha   <- if(intcens == 0L) matrix(0e0, k, 2) else NULL
  accum_row  <- NULL; accum_col <- NULL; accum_ai <- NULL; accum_ne <- 0L
  accum_hb   <- matrix(0e0, p, p)
  accum_hab  <- matrix(0e0, k, p)

  for(m in 1 : nAGQ) {
    gamma_node_c  <- gamma_c + step_c * gh$nodes[m]
    gamma_node    <- gamma_node_c[cluster]
    node_off  <- base_lp - xbeta + gamma_node        # ormll adds x%*%beta itself
    node_wt   <- wt * rho[cluster, m]
    a_arg     <- if(intcens == 0L) matrix(0e0, k, 2) else matrix(0e0, 0, 2)
    wc <- .Fortran(F_ormll, n, k, p, x, ia, ia2, sgn, ib, nb,
                   node_off, node_wt, zeromat,
                   link=link, alpha, beta, logL=numeric(1),
                   grad=numeric(k + p), lpe=numeric(n),
                   a=a_arg, b=matrix(0e0, p, p), ab=matrix(0e0, k, p),
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

  if(p > 0) {
    accum_grad[(k + 1) : (k + p)] <- accum_grad[(k + 1) : (k + p)] -
                                     as.vector(penmat %*% beta)
    if(penhess > 0) accum_hb <- accum_hb - penmat
  }

  list(fail=FALSE, loglik=total_ll, rho=rho, gh=gh, step_c=step_c,
       gamma_c=gamma_c, kappa_c=kappa_c,
       grad=accum_grad, ha=accum_ha, hb=accum_hb, hab=accum_hab,
       row=accum_row, col=accum_col, ai=accum_ai, ne=accum_ne)
}

## Cheap (ormeta-only, no ormll calls) AGQ log-likelihood at a trial
## sigma, re-finding the per-cluster mode at that sigma. Used only for
## the 1-D profile search over sigma -- never touches the O(k)/O(p^2)
## machinery.
## Claude Sonnet 5 2026-08-30          entire function
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

## Sparse missing-information correction (Louis 1982), computed via
## direct per-cluster analytic accumulation rather than global finite
## differences -- see the design discussion accompanying this file.
## Never materializes anything O(k^2): the alpha-alpha block is built
## as a sparse triplet accumulation, one small dense outer product per
## cluster (size = that cluster's own active alpha-index count, plus
## p for beta), keyed by that cluster's own alpha indices -- however
## large k is. Total nonzero count is bounded by roughly
## nc * (2*max cluster size)^2, independent of k, for realistic
## (bounded cluster size) configurations.
##
## Requires ormeta's s1/s2 outputs (the per-observation score credited
## separately to alpha(ia(i)) and alpha(ia2(i)) -- g = s1 + s2 is what
## the rest of this file uses, but the missing-information correction
## needs the two pieces kept apart, since which alpha index each piece
## credits is exactly what determines the sparsity pattern).
##
## Returns the (k+p) x (k+p) sparse missing-information matrix only;
## the caller combines it with the (already sparse, tri-band for
## alpha) complete-data information from agqStep's own ha/hb/hab.
## Claude Sonnet 5 2026-08-30          entire function
## Claude Sonnet 5 2026-08-30
## Extended to also produce the log(sigma) (tau) row/column of the
## information matrix. Writing tau = log(sigma), the prior contributes
## d/dtau[log phi(u;0,sigma)] = -1 + u^2/sigma^2 to the complete-data
## score, and this does NOT depend on theta=(alpha,beta) at all (the
## data likelihood given u doesn't involve sigma) -- so the complete-
## data alpha/beta-tau cross information is EXACTLY zero, and all
## alpha/beta-tau structure in the final (corrected) information
## matrix comes purely from the missing-information term below, via
## one more column appended to each cluster's own score vector. The
## complete-data tau-tau diagonal entry (from d2/dtau2 of the prior
## term, = -2u^2/sigma^2, expectation taken via the same AGQ weights)
## is accumulated in the same per-node loop already computing
## everything else needed here, and returned separately for the
## caller to fold into the complete-data block directly (it is not,
## and cannot be, part of the sparse per-cluster alpha/beta
## accumulation, since it is a single scalar spanning every cluster).
sparseMissingInfo <- function(alpha, beta, sigma, gamma_c, kappa_c, base_lp,
                              x, wt, ia, ia2, sgn, ib, nb,
                              cluster, nc, link, k, p, n, nAGQ) {
  gh     <- gauss_hermite_quad(nAGQ)
  step_c <- sqrt(2 / kappa_c)
  ptau   <- p + 1L   # p fixed-effect columns plus one for log(sigma)

  obs_by_cluster <- split(seq_len(n), cluster)

  logd_nodes <- matrix(0e0, n, nAGQ)
  s1_nodes   <- matrix(0e0, n, nAGQ)
  s2_nodes   <- matrix(0e0, n, nAGQ)
  g_nodes    <- matrix(0e0, n, nAGQ)
  gamma_node_mat_full <- matrix(0e0, nc, nAGQ)   # per-cluster node values, all nodes

  for(m in 1 : nAGQ) {
    gamma_node_c <- gamma_c + step_c * gh$nodes[m]
    gamma_node_mat_full[, m] <- gamma_node_c
    gamma_node   <- gamma_node_c[cluster]
    lp       <- base_lp + gamma_node
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

  ## Per-cluster, per-node AGQ log-lik and posterior node weights
  clusterloglik <- matrix(0e0, nc, nAGQ)
  for(j in 1 : nc) {
    obs <- obs_by_cluster[[j]]
    clusterloglik[j, ] <- colSums(logd_nodes[obs, , drop=FALSE])
  }
  lt   <- sweep(clusterloglik, 2, log(gh$weights) + gh$nodes^2, "+") +
          dnorm(gamma_node_mat_full, mean=0, sd=sigma, log=TRUE)
  mrow <- apply(lt, 1, max)
  rho  <- exp(lt - mrow); rho <- rho / rowSums(rho)   # nc x nAGQ

  ## Complete-data tau-tau diagonal: E[-d2/dtau2 of the prior term |
  ## data] = 2/sigma^2 * E[u^2 | data], AGQ-weighted, summed over
  ## clusters. Pure R, no Fortran needed -- purely a function of
  ## gamma_node_mat_full/rho/sigma already in hand.
  complete_tau_info <- sum(rho * gamma_node_mat_full^2) * 2 / sigma^2

  ## Per-cluster small dense outer products, accumulated into a
  ## global sparse triplet list. Dimension now na + p + 1 (the extra
  ## 1 for tau); tau's own per-node score contribution is
  ## -1 + gamma_node^2/sigma^2 -- the constant -1 contributes nothing to
  ## the AGQ-weighted covariance (already subtracted via Sbarj below)
  ## but is included for clarity/direct correspondence to the math.
  tI <- integer(0); tJ <- integer(0); tX <- numeric(0)

  for(j in 1 : nc) {
    obs    <- obs_by_cluster[[j]]
    iao    <- ia[obs]; ia2o <- ia2[obs]
    active <- sort(unique(c(iao, ia2o[ia2o > 0])))
    na     <- length(active)
    dimj   <- na + ptau   # always >= 1 now (tau is always present)

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
      Sj[m, dimj] <- -1 + gamma_node_mat_full[j, m]^2 / sigma^2   # tau column
    }
    wj    <- rho[j, ]
    Sbarj <- as.vector(wj %*% Sj)
    Dj    <- sweep(Sj, 2, Sbarj)
    Vj    <- crossprod(Dj, wj * Dj)   # AGQ-weighted covariance, dimj x dimj

    idxglobal <- c(active, if(p > 0) (k + 1) : (k + p) else integer(0), k + p + 1)
    tI <- c(tI, rep(idxglobal, times = dimj))
    tJ <- c(tJ, rep(idxglobal, each  = dimj))
    tX <- c(tX, as.vector(Vj))
  }

  missing_info <- if(length(tI) > 0)
    Matrix::sparseMatrix(i=tI, j=tJ, x=tX, dims=c(k + p + 1, k + p + 1)) else
    Matrix::Matrix(0, k + p + 1, k + p + 1, sparse=TRUE)

  list(missing_info = missing_info, complete_tau_info = complete_tau_info)
}

## Assembles the corrected, sparse (k+p) x (k+p) information matrix:
## the existing (already sparse, tri-band for alpha) complete-data
## information from agqStep's ha/hb/hab, minus the missing-information
## correction above (Louis 1982). The optimization loop's own use of
## ha/hb/hab (picking a Newton direction every iteration) is entirely
## untouched -- this assembly happens once, after convergence.
## Claude Sonnet 5 2026-08-30          entire function
## Claude Sonnet 5 2026-08-30          restructured return value
## Returns info in EXACTLY the format ormfit itself uses for its
## interval-censoring case: info$a is a (row, col, a) triplet list
## (general sparse, upper-triangle-only per Matrix::sparseMatrix's
## symmetric=TRUE convention) rather than the k x 2 tri-band format
## used when no correction is needed -- appropriate here because
## integrating out a shared random intercept can couple non-adjacent
## thresholds regardless of whether the underlying data has interval
## censoring (which orm.rfit doesn't yet support anyway). info$b/ab
## are plain dense matrices, exactly as ormfit already returns.
## This makes infoMxop(f$info.matrix, invert=TRUE) work unchanged --
## no separate extraction path needed for orm.rfit's output.
sparseInfoMatrix <- function(alpha, beta, sigma, gamma_c, kappa_c, base_lp,
                             x, wt, ia, ia2, sgn, ib, nb, cluster, nc,
                             link, k, p, n, ha, hb, hab, row, col, ai, ne,
                             intcens, nAGQ, iname, xname) {
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

  ## Claude Sonnet 5 2026-08-30
  ## Extend hb (p x p) to (p+1) x (p+1) and hab (k x p) to k x (p+1)
  ## for log(sigma) (tau): complete-data alpha/beta-tau cross terms
  ## are exactly zero (see sparseMissingInfo's header comment), so
  ## the extra row/column is zero except the tau-tau diagonal entry,
  ## which sparseMissingInfo computed alongside everything else it
  ## already needed (complete_tau_info, on the raw Hessian scale, so
  ## negated below along with everything else for the info-matrix
  ## sign convention).
  ptau      <- p + 1L
  hb_ext    <- matrix(0e0, ptau, ptau)
  if(p > 0) hb_ext[1 : p, 1 : p] <- hb
  hab_ext   <- matrix(0e0, k, ptau)
  if(p > 0) hab_ext[, 1 : p] <- hab

  mi <- sparseMissingInfo(alpha, beta, sigma, gamma_c, kappa_c, base_lp,
                          x, wt, ia, ia2, sgn, ib, nb,
                          cluster, nc, link, k, p, n, nAGQ)
  hb_ext[ptau, ptau] <- -mi$complete_tau_info   # raw Hessian scale, negated below with the rest

  complete_data_info <-
    rbind(cbind(-Ha_sparse,                -Matrix::Matrix(hab_ext, sparse=TRUE)),
          cbind(-Matrix::Matrix(t(hab_ext), sparse=TRUE), -Matrix::Matrix(hb_ext, sparse=TRUE)))

  V <- complete_data_info - mi$missing_info

  ## Split into the a (sparse triplet, upper-triangle only) / b / ab
  ## pieces ormfit's own intcens=1 output uses. b is now (p+1)x(p+1)
  ## and ab is k x (p+1), with the last column/row being log(sigma).
  Va <- V[1 : k, 1 : k, drop=FALSE]
  Vsumm <- Matrix::summary(methods::as(Matrix::drop0(Va), "TsparseMatrix"))
  keep  <- Vsumm$i <= Vsumm$j
  a_list <- list(row = Vsumm$i[keep], col = Vsumm$j[keep], a = Vsumm$x[keep])

  b_mat  <- as.matrix(V[(k + 1) : (k + ptau), (k + 1) : (k + ptau), drop=FALSE])
  ab_mat <- as.matrix(V[1 : k, (k + 1) : (k + ptau), drop=FALSE])

  list(a = a_list, b = b_mat, ab = ab_mat,
      iname = iname, xname = c(xname, "log(sigma)"))
}

ormrfit <- function(x, y, y2, k, intcens=0L, cluster, nc, initial, sigma.init=1.0,
                    offset=rep(0., n), wt=rep(1., n), penmat=matrix(0., p, p),
                    maxit=30L, maxit.outer=100L, maxit.mode=30L,
                    objtol=5e-4, gradtol=1e-3, paramtol=1e10,
                    tolsolve=.Machine$double.eps, minstepsize=1e-2, trace=FALSE,
                    link, iname, xname,
                    nAGQ=7L, nAGQ.grid=c(7L, 11L, 15L, 21L, 31L, 45L, 63L),
                    nAGQ.tol=1e-5) {

  # Claude Sonnet 5 2026-08-30          entire function
  n <- length(y)
  p <- length(initial) - k

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
  gamma_c     <- rep(0e0, nc)
  xbeta   <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)
  base_lp <- offset + xbeta

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
    theta_full <- c(theta, log(sigma))

    ## Claude Sonnet 5 2026-09-04          entire onestep() function
    ## One full outer "fixed-point" iteration -- per-cluster
    ## mode-finding, the AGQ-reweighted Newton step for (alpha,beta)
    ## with its own step-halving, and the sigma profile search --
    ## factored out unchanged from the previous plain loop so it can
    ## be called repeatedly per SQUAREM cycle below. grad/delta in its
    ## return value are evaluated AT theta_full going IN (matching the
    ## pre-existing semantics exactly: gradient before the step, size
    ## of that step), not at the result.
    onestep <- function(theta_full, gamma_c, oi, role) {
      alpha   <- theta_full[1 : k]
      beta    <- if(p > 0) theta_full[(k + 1) : (k + p)] else numeric(0)
      sigma   <- exp(theta_full[k + p + 1])
      theta   <- theta_full[1 : (k + p)]
      xbeta   <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)
      base_lp <- offset + xbeta

      mf <- clusterModeFind(alpha, beta, sigma, gamma_c, base_lp,
                            wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                            cluster, nc, link, k, n, maxit=maxit.mode)
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

      ag <- agqStep(alpha, beta, sigma, gamma_c, mf$kappa_c, cur.nAGQ,
                   base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                   cluster, nc, link, k, p, n, penmat, penhess=1L,
                   intcens=intcens, score=TRUE)
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
      hess  <- infoMxop(list(a=a_for_hess, b=ag$hb, ab=ag$hab))
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
      lambda <- 0
      repeat {
        hess_try <- if(lambda == 0) hess else
                      hess + lambda * Matrix::Diagonal(x=Matrix::diag(hess))
        delta <- try(Matrix::solve(hess_try, ag$grad, tol=tolsolve), silent=TRUE)
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
      surrogate <- function(theta_t) {
        alpha_t <- theta_t[1 : k]
        if(k > 1 && any(diff(alpha_t) >= 0)) return(Inf)
        beta_t    <- if(p > 0) theta_t[(k + 1) : (k + p)] else numeric(0)
        xbeta_t   <- if(p > 0) as.vector(x %*% beta_t) else rep(0e0, n)
        lp_base_t <- offset + xbeta_t
        tot <- 0
        for(mm in 1 : cur.nAGQ) {
          gamma_node_c <- gamma_c + ag$step_c * ag$gh$nodes[mm]
          gamma_node   <- gamma_node_c[cluster]
          lp_t     <- lp_base_t + gamma_node
          wtt <- .Fortran(F_ormeta, n, k, link, alpha_t, lp1=lp_t, lp2=lp_t, wt,
                          idx$ia, idx$ia2, idx$sgn, ib, nb,
                          logd=numeric(n), g=numeric(n), h=numeric(n),
                          s1=numeric(n), s2=numeric(n),
                          debug=0L, salloc=integer(1))
          if(wtt$salloc != 0) return(Inf)
          tot <- tot + sum(ag$rho[cluster, mm] * wtt$logd)
        }
        -2 * tot
      }

      objf      <- surrogate(theta)   # current theta, same fixed rho
      step_size <- 1.0
      repeat {
        theta_new <- theta - step_size * delta
        objfnew   <- surrogate(theta_new)
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

      alpha_new <- theta_new[1 : k]
      beta_new  <- if(p > 0) theta_new[(k + 1) : (k + p)] else numeric(0)
      xbeta_new <- if(p > 0) as.vector(x %*% beta_new) else rep(0e0, n)
      base_lp_new <- offset + xbeta_new

      ## sigma update: 1-D profile search, ormeta-only (cheap)
      br   <- log(sigma) + c(-2.5, 2.5)
      sopt <- optimize(sigmaObjective, interval=br, maximum=TRUE,
                       alpha=alpha_new, beta=beta_new, gamma_c_init=gamma_c, base_lp=base_lp_new,
                       wt=wt, ia=idx$ia, ia2=idx$ia2, sgn=idx$sgn, ib=ib, nb=nb,
                       cluster=cluster, nc=nc, link=link, k=k, n=n,
                       nAGQ=cur.nAGQ, maxit.mode=maxit.mode, tol=1e-4)
      sigma_new <- exp(sopt$maximum)
      sigma.change <- abs(sigma_new - sigma) / max(sigma, 1e-6)
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
            ' max|grad|:', m(ag$grad), ' max|delta|:', m(delta), ']\n')

      list(fail=FALSE, theta_full=c(theta_new, log(sigma_new)),
          gamma_c=gamma_c, grad=ag$grad, delta=delta,
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
      if(trace)
        cat('nAGQ:', cur.nAGQ, ' outer iter:', outer.iter,
            ' sigma:', format(exp(theta_full[k + p + 1]), nsmall=4),
            ' max|grad|:', m(accepted$grad), ' max|delta|:', m(accepted$delta), '\n')

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
        msg <- paste('Reached', maxit.outer, 'outer iterations without convergence at nAGQ=',
                    cur.nAGQ, '\nsigma:', exp(theta_full[k + p + 1]),
                    ' Max |gradient|:', m(accepted$grad),
                    ' Max |change in parameters|:', m(accepted$delta))
        message(msg)
        return(list(fail=TRUE))
      }
    }

    total.oi <- total.oi + outer.iter

    alpha   <- theta_full[1 : k]
    beta    <- if(p > 0) theta_full[(k + 1) : (k + p)] else numeric(0)
    sigma   <- exp(theta_full[k + p + 1])
    theta   <- theta_full[1 : (k + p)]
    xbeta   <- if(p > 0) as.vector(x %*% beta) else rep(0e0, n)
    base_lp <- offset + xbeta

    ## Escalation check: does the next larger nAGQ change the
    ## log-likelihood by more than nAGQ.tol? (cheap, ormeta-only,
    ## reusing the mode/theta just converged to as a warm start)
    mf.cur <- clusterModeFind(alpha, beta, sigma, gamma_c, base_lp,
                              wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                              cluster, nc, link, k, n, maxit=maxit.mode)
    ll.cur <- agqStep(alpha, beta, sigma, mf.cur$gamma_c, mf.cur$kappa_c, cur.nAGQ,
                      base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                      cluster, nc, link, k, p, n, penmat, score=FALSE)$loglik

    if(grid.pos == length(nAGQ.grid)) {
      if(trace) message('nAGQ reached grid ceiling (', cur.nAGQ,
                        ') without full convergence -- results may still be changing')
      break
    }
    next.nAGQ <- nAGQ.grid[grid.pos + 1]
    ll.next   <- agqStep(alpha, beta, sigma, mf.cur$gamma_c, mf.cur$kappa_c, next.nAGQ,
                        base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                        cluster, nc, link, k, p, n, penmat, score=FALSE)$loglik
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
  ## the returned information matrix
  mf.final <- clusterModeFind(alpha, beta, sigma, gamma_c, base_lp,
                              wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                              cluster, nc, link, k, n, maxit=maxit.mode)
  ag.final <- agqStep(alpha, beta, sigma, mf.final$gamma_c, mf.final$kappa_c, cur.nAGQ,
                      base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                      cluster, nc, link, k, p, n, penmat,
                      intcens=intcens, score=TRUE)

  info <- sparseInfoMatrix(alpha, beta, sigma, mf.final$gamma_c, mf.final$kappa_c,
                           base_lp, x, wt, idx$ia, idx$ia2, idx$sgn, ib, nb,
                           cluster, nc, link, k, p, n,
                           ag.final$ha, ag.final$hb, ag.final$hab,
                           ag.final$row, ag.final$col, ag.final$ai, ag.final$ne,
                           intcens, cur.nAGQ, iname, xname)

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
  list(coef    = theta,
      sigma    = sigma,
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
