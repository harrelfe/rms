predict.lrm <- function(object, ..., 
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		"terms", "cterms", "ccterms", "adjto", "adjto.data.frame",
          "model.frame"),
		se.fit=FALSE, codes=FALSE)
{
  type <- match.arg(type)
  if(type %nin% c("fitted","fitted.ind", "mean"))
    return(predictrms(object,...,type=type, se.fit=se.fit))

  xb <- predictrms(object, ..., type="lp", se.fit=FALSE)
  rnam <- names(xb)
  ns <- object$non.slopes
  cnam <- names(object$coef[1:ns])
  cumprob <- if(length(object$famfunctions)) eval(object$famfunctions[1]) else plogis
  if(se.fit)
    warning('se.fit not supported with type="fitted" or type="mean"')
  if(ns == 1 & type == "mean")
    stop('type="mean" makes no sense with a binary response')
  if(ns == 1) return(cumprob(xb))
  intcept <- object$coef[1:ns]
  interceptRef <- object$interceptRef
  if(!length(interceptRef)) interceptRef <- 1
  xb <- xb - intcept[interceptRef]
  xb <- sapply(intcept, "+", xb)
  P <- cumprob(xb)
  nam <- names(object$freq)
  if(is.matrix(P)) dimnames(P) <- list(rnam, cnam)
  else names(P) <- names(object$coef[1:ns])
  if(type=="fitted") return(P)

  ##type="mean" or "fitted.ind"
  vals <- object$yunique
  P   <- matrix(P, ncol=ns)
  Peq <- cbind(1, P) - cbind(P, 0)
  if(type == "fitted.ind") {
    ynam <- as.character(attr(object$terms, "formula")[2])
    ynam <- paste(ynam, "=", vals, sep="")
    dimnames(Peq) <- list(rnam, ynam)
    return(drop(Peq))
  }
  
  ##type="mean"
  if(codes) vals <- 1:length(object$freq)
  else {
    vals <- as.numeric(vals)
    if(any(is.na(vals)))
      stop('values of response levels must be numeric for type="mean" and codes=F')
  }
  m <- drop(Peq %*% vals)
  names(m) <- rnam
  m
}

predict.orm <- function(object, ...,
		type=c("lp","fitted","fitted.ind","mean","x","data.frame",
		"terms", "cterms", "ccterms", "adjto", "adjto.data.frame",
          "model.frame"),
		se.fit=FALSE, codes=FALSE)
{
  type <- match.arg(type)
  # Claude Sonnet 5 2026-11-20          8 lines
  # A yde-fitted orm() object (qyde > 0) needs the Y-dependent (tau)
  # contribution included -- predict.lrm()/predictrms() know nothing
  # about it (no modifications to lrm() or its own methods; see
  # orm-yde-design.md section b.4). Route EVERY type through
  # ormPredictYde() when yde is active, not just the types it happens to
  # implement -- letting an unimplemented type fall through to
  # predict.lrm()/predictrms() would silently build a linear predictor
  # missing the Y-dependent contribution instead of erroring (confirmed
  # this was happening for type='fitted'/'mean' before this fix).
  # ormPredictYde() itself decides what's implemented vs. stops.
  if(length(object$qyde) && object$qyde > 0)
    return(ormPredictYde(object, ..., type=type, se.fit=se.fit, codes=codes))
  predict.lrm(object, ..., type=type, se.fit=se.fit, codes=codes)
}

##' Y-Dependent-Effect-Aware Prediction for an `orm()` Fit Using `yde`/`ytrans`
##'
##' Internal helper for [predict.orm()]: implements `type='lp'` and
##' `type='x'` for a `yde`-fitted `orm()` object, adding the Y-dependent
##' (`tau`) contribution that [predict.lrm()]/[predictrms()] know
##' nothing about. Reuses [predictrms()]'s existing `second=TRUE`
##' mechanism (built for [rmsb::blrm()]'s own `ppo`, reading
##' `object$zDesign`/`object$zassign`/`object$zsformula`) to build the
##' `yde`-only design matrix from `newdata`, rather than a second,
##' hand-rolled `Design()`/`model.matrix()` implementation -- the
##' analog, for an MLE fit, of [bormPredictYde()] in `rmsb`. `kint`/
##' `ycut` indexing follows the same convention already used by
##' [Predict()] and [rmsb::predict.blrm()]: `kint <- which(yunique ==
##' ycut) - 1` lands on the same 1-based index that both `coef(object)`'s
##' alpha slice uses (alpha index `i` <-> `yunique[i + 1]`, confirmed
##' from how `iname`/`ylev` are built in `orm.fit`). `object$ytransmat`
##' itself has `k+1` rows, not `k` (Claude Sonnet 5 2026-11-23, matching
##' `rmsb::borm()`'s own convention): row 1 is the baseline Y level
##' (`yunique[1]`, with no intercept of its own), so alpha index `i`'s
##' own row is `ytransmat[i + 1, ]`, not `ytransmat[i, ]`.
##' @param object an `orm()` fit object with `qyde > 0`
##' @param type `'lp'` or `'x'`; nothing else is implemented yet
##' @param kint intercept index (1-based); see Details
##' @param ycut alternative to `kint`: a Y value to evaluate the
##'   Y-dependent effect at
##' @author Frank Harrell
##' @noRd
ormPredictYde <- function(object, ..., type, kint=NULL, ycut=NULL, se.fit=FALSE,
                           Zmatrix=TRUE, zcppo=TRUE, fun=NULL, funint=TRUE,
                           codes=FALSE, posterior.summary=NULL, cint=NULL) {
  # Claude Sonnet 5 2026-11-20          4 lines
  # 'fitted'/'fitted.ind'/'mean' are handled by a separate helper,
  # ormYdeFittedMean() (below) -- a different enough computation (all k
  # cutpoints at once, or streamed over them) that threading it through
  # this function's existing single-ycut lp/x logic would be more
  # confusing than two focused functions sharing the same setup.
  implemented <- c('lp', 'x', 'fitted', 'fitted.ind', 'mean')
  if(type %nin% implemented)
    stop("type='", type, "' is not yet implemented for an orm() fit using yde/ytrans; ",
         "only ", paste(sQuote(implemented), collapse=', '),
         " are currently supported for such fits")
  if(se.fit && type == 'x')
    warning('se.fit does not apply to type="x"')
  # Claude Sonnet 5 2026-11-20
  # type='mean' now has a real delta-method se.fit (see
  # ormYdeFittedMean()) -- only 'fitted'/'fitted.ind' still ignore it,
  # matching predict.lrm()'s own (non-yde) convention for those two.
  if(se.fit && type %in% c('fitted', 'fitted.ind'))
    warning('se.fit not supported with type="fitted" or "fitted.ind"')
  if(length(fun))
    stop('fun= is not yet implemented for an orm() fit using yde/ytrans')
  # Claude Sonnet 5 2026-11-17
  # funint/codes/posterior.summary/cint are consumed here, as explicit
  # formals, purely so they do NOT fall through to ... and get forwarded
  # to predictrms() below -- predictrms()'s own model.frame() call
  # treats any unrecognized named ... argument as an extra variable to
  # evaluate in newdata, so passing fun=NULL/funint=TRUE etc. straight
  # through (as a first attempt did) crashes with "object 'fun' not
  # found", not a silent no-op. None of these apply to an MLE fit.

  X <- predictrms(object, ..., type='x')
  # Claude Sonnet 5 2026-11-17          2 lines
  # Zmatrix=FALSE: matching predict.blrm()/bormPredictYde()'s own early
  # exit, skip the (second predictrms() call plus) yde design-matrix
  # construction below entirely. zcppo accepted for interface
  # compatibility only, same as bormPredictYde() -- see its own
  # documented reason (tau has no meaningful "unweighted" form).
  if(type == 'x' && ! Zmatrix) return(X)
  Z <- predictrms(object, ..., type='x', second=TRUE)
  n <- NROW(X)

  if(type %in% c('fitted', 'fitted.ind', 'mean'))
    return(ormYdeFittedMean(object, X=X, Z=Z, type=type, codes=codes,
                             se.fit = se.fit && type == 'mean'))

  ylevels   <- object$yunique
  ytransmat <- object$ytransmat
  qyde      <- object$qyde
  myde      <- ncol(ytransmat)

  iref <- object$interceptRef
  if(! length(kint)) kint <- iref
  ycutgiven <- length(ycut) > 0
  if(! ycutgiven) ycut <- ylevels[kint + 1]
  if(ycutgiven) {
    if(type == 'lp' && length(ycut) > 1) stop('ycut may only be a scalar for type="lp"')
    # Claude Sonnet 5 2026-11-22          entire block
    # Fixed bug: for character/non-numeric ylevels this used to do
    # max((1:length(ylevels))[ylevels <= ycut[1]]) - 1, an ALPHABETICAL
    # string comparison, not a lookup by the levels' actual ordinal
    # position. Confirmed directly (reported by Frank from a real
    # fit): with ylevels = c("Home", "In Hospital/Facility",
    # "Vent/ARDS", "Dead"), "Dead" sorts alphabetically before all the
    # others (D < H < I < V), so it is always included in the <= set
    # and always wins the max() -- giving the SAME (wrong) index for
    # ycut="In Hospital/Facility", ycut="Dead", and even a nonexistent
    # ycut="junk", which is why none of these behaved differently and
    # why an invalid ycut never errored. Replaced with match(), an
    # exact lookup by position in ylevels (which IS already in the
    # correct modeled order) -- the same exact-match semantics the
    # numeric branch already used (which(ylevels==ycut[1])), just
    # written so it also works for character/factor ylevels, with no
    # separate branch needed. An unmatched ycut now stops with an
    # informative error instead of silently resolving to some other
    # level's index.
    kint <- match(ycut[1], ylevels) - 1
    if(is.na(kint))
      stop('ycut=', ycut[1], ' does not match any modeled Y level; valid levels are: ',
           paste(ylevels, collapse=', '))
  }

  if(length(ycut) %nin% c(1, n))
    stop('ycut must be of length 1 or the number of requested predictions')
  ycut   <- rep(ycut, length = n)
  rowidx <- match(ycut, ylevels) - 1L
  if(any(is.na(rowidx)))
    stop('ycut has value(s) not matching any modeled Y level: ',
         paste(unique(ycut[is.na(rowidx)]), collapse=', '),
         '; valid levels are: ', paste(ylevels, collapse=', '))
  if(any(rowidx < 1L))
    stop('ycut below the lowest modeled cutpoint has no defined Y-dependent effect')

  expandZ <- function(Zrow, ytransrow) as.vector(outer(Zrow, ytransrow))
  Zexp <- matrix(0., n, qyde * myde)
  # Claude Sonnet 5 2026-11-23
  # rowidx is a 1-based alpha index (1..k); ytransmat now has k+1 rows
  # (row 1 is the baseline Y level, with no intercept of its own), so
  # row (rowidx + 1) is the one that corresponds to a given alpha
  # index -- see orm.fit.r's own identical shift on idx$ia.
  for(i in 1 : n) Zexp[i, ] <- expandZ(Z[i, ], ytransmat[rowidx[i] + 1L, ])
  taunam <- utils::tail(names(coef(object)), qyde * myde)
  colnames(Zexp) <- taunam

  if(type == 'x') return(cbind(X, Zexp))

  ## type == 'lp'
  k    <- object$non.slopes
  p    <- ncol(X)
  qm   <- qyde * myde
  cof  <- coef(object)
  alpha <- cof[rowidx]
  beta  <- cof[(k + 1) : (k + p)]
  tau   <- cof[(k + p + 1) : (k + p + qm)]
  lp <- alpha + (if(p > 0) as.vector(X %*% beta) else 0.) + as.vector(Zexp %*% tau)
  names(lp) <- rownames(X)

  # Claude Sonnet 5 2026-11-17          entire block
  # Delta-method SE: v_i = [one-hot at alpha index rowidx[i] (length
  # k), X[i,], Zexp[i,]], se_i = sqrt(v_i' Vsub v_i), where Vsub is the
  # covariance submatrix for the *distinct* alpha indices actually
  # touched plus all of beta+tau -- requested from infoMxop() by name
  # (its own efficient, non-full-inversion path for an arbitrary index
  # set) rather than inverting the full (possibly huge, k large) info
  # matrix. Needed because Predict()'s default conf.int=0.95 means
  # Pred()'s caller unconditionally reads xx$lower/xx$upper -- without
  # this, a yde-fitted orm() object crashed there ("$ operator is
  # invalid for atomic vectors") on every default-options Predict()
  # call, not just when the user asked for se.fit/conf.int explicitly.
  if(! se.fit && ! length(cint)) return(lp)
  uniq_a <- sort(unique(rowidx))
  idxset <- c(uniq_a, (k + 1) : (k + p + qm))
  Vsub   <- infoMxop(object$info.matrix, i = idxset)
  napos  <- length(uniq_a)
  apos   <- match(rowidx, uniq_a)
  bpos   <- (napos + 1) : (napos + p + qm)
  se <- numeric(n)
  for(i in 1 : n) {
    v <- numeric(napos + p + qm)
    v[apos[i]] <- 1.
    if(p > 0) v[bpos[1 : p]] <- X[i, ]
    v[bpos[(p + 1) : (p + qm)]] <- Zexp[i, ]
    se[i] <- sqrt(drop(v %*% Vsub %*% v))
  }
  if(! length(cint)) return(list(linear.predictors = lp, se.fit = se))
  zcrit <- qnorm((1. + cint) / 2.)
  list(linear.predictors = lp, se.fit = se,
       lower = lp - zcrit * se, upper = lp + zcrit * se)
}

##' Fitted/Mean Values for a Y-Dependent-Effect `orm()` Fit
##'
##' Internal helper for [ormPredictYde()]: implements `type='fitted'`,
##' `'fitted.ind'`, and `'mean'` for a `yde`-fitted `orm()` object.  All
##' three need the Y-dependent contribution at *every* cutpoint `j`
##' (1..k), unlike [ormPredictYde()]'s own `'lp'`/`'x'` branches, which
##' only ever need one (`ycut`/`kint`). Exploits the same bilinear
##' structure for that: writing `tau` as a `qyde x myde` matrix `Tau`,
##' the Y-dependent contribution at cutpoint `j` is
##' `Z %*% Tau %*% ytransmat[j + 1, ]` (Claude Sonnet 5 2026-11-23:
##' `+ 1` because `ytransmat`'s own row 1 is the baseline Y level, not
##' cutpoint 1 -- see [ormPredictYde()]'s own doc comment), which is
##' bilinear in `i` (observation) and `j` (cutpoint) and so factors
##' into two pieces computed once, independent of `k`: `Xbeta <- X %*%
##' beta` and `W <- Z %*% Tau` (`n x myde`). `'fitted'`/`'fitted.ind'`
##' then build the full `n x k` matrix of cutpoint linear predictors in
##' one matrix product (`W %*% t(ytransmat[-1, ])`), exactly as fast
##' asymptotically as just
##' materializing the `n x k` output requires (confirmed by complexity
##' analysis, not merely assumed -- see `orm-yde-implementation.md`), and
##' difference the resulting cumulative probabilities exactly as
##' [predict.lrm()] already does for an ordinary (non-`yde`) fit.
##' `'mean'` instead uses an algebraic identity -- verified directly by
##' expanding and telescoping [predict.lrm()]'s own
##' `(cbind(1,P)-cbind(P,0)) %*% values` (a general property of
##' category-probabilities-from-cumulative-probabilities, not special to
##' `yde`): `E[Y] = values[1] + sum_j (values[j+1]-values[j]) * Pcum[,j]`
##' -- to accumulate the mean via a loop over the `k` cutpoints that
##' only ever holds an `n`-length vector, rather than ever materializing
##' the full `n x k` matrix. That trades nothing in speed (`cumprob`
##' still genuinely runs `n*k` times either way) for `O(n)` instead of
##' `O(n*k)` peak memory -- the lever that matters for `orm()`'s
##' large-`k` regime. This memory saving is specific to a *scalar*
##' summary like the mean, where the `k` dimension collapses
##' algebraically before anything needs to be stored; it does not apply
##' to `'fitted'`/`'fitted.ind'`, which are asked to return the full
##' per-category table.
##' @param object a `yde`-fitted `orm()` object
##' @param X the `n x p` main-effect design matrix, as returned by
##'   `predictrms(object, ..., type='x')`
##' @param Z the `n x qyde` `yde` design matrix, as returned by
##'   `predictrms(object, ..., type='x', second=TRUE)`
##' `type='mean'`'s `se.fit=TRUE` path is a genuinely different, separate
##' code path from its point-estimate-only one (Frank's explicit request:
##' keep the fastest, `O(n)` computation for the common case where SE
##' isn't wanted, even though that means two code paths instead of one).
##' It needs the full `n x k` matrix regardless, since the delta-method
##' gradient has a nonzero entry for every alpha (every cutpoint
##' contributes to the mean, unlike `'lp'`/`'fitted.ind'`, where
##' [ormPredictYde()]'s own `se.fit` code gets away with a small "touched
##' indices" subset) -- not a new cost `yde` introduces: [Mean.lrm()]'s
##' existing (non-`yde`) `conf.int` code already pays this same `O(n*k)`
##' cost via the same `infoMxop(..., B=)` mechanism reused here, just
##' with `qm` extra `tau` columns appended to `dmean.dtheta`, computed via
##' the same bilinear identity as the point estimate
##' (`d(mean)/d(tau) = Z \%o\% (dmean.dalpha \%*\% ytransmat[-1, ])`,
##' verified by direct differentiation of the telescoping identity
##' above; `ytransmat`'s own row 1, the baseline Y level, is dropped
##' for the same reason as elsewhere in this function).
##' @param object a `yde`-fitted `orm()` object
##' @param X the `n x p` main-effect design matrix, as returned by
##'   `predictrms(object, ..., type='x')`
##' @param Z the `n x qyde` `yde` design matrix, as returned by
##'   `predictrms(object, ..., type='x', second=TRUE)`
##' @param type `'fitted'`, `'fitted.ind'`, or `'mean'`
##' @param codes see [predict.lrm()]'s own `codes` argument (only used
##'   by `type='mean'`)
##' @param se.fit only used by `type='mean'`; delta-method standard
##'   error of the mean (see Details)
##' @author Frank Harrell
##' @noRd
ormYdeFittedMean <- function(object, X, Z, type, codes=FALSE, se.fit=FALSE) {
  ylevels   <- object$yunique
  ytransmat <- object$ytransmat
  qyde      <- object$qyde
  myde      <- ncol(ytransmat)
  k         <- object$non.slopes
  p         <- ncol(X)
  qm        <- qyde * myde
  cof       <- coef(object)
  alpha     <- cof[1 : k]
  beta      <- cof[(k + 1) : (k + p)]
  tau       <- cof[(k + p + 1) : (k + p + qm)]
  # Claude Sonnet 5 2026-11-20
  # Tau's flattening order matches ormPredictYde()'s own expandZ()/
  # as.vector(outer(Zrow, ytransrow)) convention (column-major by
  # ytrans column), so matrix(tau, qyde, myde) (default byrow=FALSE)
  # reconstructs Tau[q, m] correctly -- confirmed, not assumed, from
  # that existing code's own comment on the flattening order.
  Tau   <- matrix(tau, nrow = qyde, ncol = myde)
  Xbeta <- if(p > 0) as.vector(X %*% beta) else 0.
  W     <- Z %*% Tau   # n x myde -- independent of k, computed once
  cumprob <- if(length(object$famfunctions)) eval(object$famfunctions[1]) else plogis
  n    <- NROW(X)
  rnam <- rownames(X)

  if(type %in% c('fitted', 'fitted.ind')) {
    # Claude Sonnet 5 2026-11-23: ytransmat's row 1 is the baseline Y
    # level (no intercept of its own) -- drop it so the remaining k
    # rows line up with alpha's own k-length indexing.
    ZTau    <- W %*% t(ytransmat[-1, , drop = FALSE])   # n x k, one matrix product
    lp_full <- sweep(ZTau, 2, alpha, '+') + Xbeta
    P       <- cumprob(lp_full)
    if(type == 'fitted') {
      dimnames(P) <- list(rnam, names(alpha))
      return(drop(P))
    }
    Peq  <- cbind(1, P) - cbind(P, 0)
    ynam <- as.character(attr(object$terms, "formula")[2])
    ynam <- paste(ynam, "=", ylevels, sep="")
    dimnames(Peq) <- list(rnam, ynam)
    return(drop(Peq))
  }

  ## type == 'mean'
  if(codes) values <- 1 : length(ylevels)
  else {
    values <- as.numeric(ylevels)
    if(any(is.na(values)))
      stop('values of response levels must be numeric for type="mean" and codes=FALSE')
  }
  dvalues <- diff(values)

  if(! se.fit) {
    # Claude Sonnet 5 2026-11-20
    # Fastest path: O(n) memory, looping over the k cutpoints instead of
    # ever materializing the full n x k matrix -- see this function's
    # own doc comment above for the telescoping identity this
    # implements. Kept as a separate code path from the se.fit=TRUE
    # block below (which needs the full matrix anyway for the
    # gradient), per Frank's explicit request to always use the fastest
    # computation when SE isn't requested, even at the cost of two code
    # paths instead of one.
    m <- numeric(n)
    # Claude Sonnet 5 2026-11-23: j+1, not j -- ytransmat row 1 is the
    # baseline Y level, so alpha[j] (1..k) lines up with ytransmat row
    # j+1 (2..k+1).
    for(j in 1 : k)
      m <- m + dvalues[j] * cumprob(alpha[j] + Xbeta + as.vector(W %*% ytransmat[j + 1L, ]))
    m <- m + values[1]
    names(m) <- rnam
    return(m)
  }

  # Claude Sonnet 5 2026-11-20          entire block
  # se.fit=TRUE needs the full n x k matrix regardless (for the
  # gradient), so the point estimate is recomputed here from it via the
  # same telescoping identity as one matrix-vector product (P %*%
  # dvalues) instead of the loop above -- no accuracy or speed lost by
  # not reusing the O(n) loop's result, since P has to be built anyway.
  ZTau    <- W %*% t(ytransmat[-1, , drop = FALSE])
  lp_full <- sweep(ZTau, 2, alpha, '+') + Xbeta
  P       <- cumprob(lp_full)
  m       <- values[1] + as.vector(P %*% dvalues)
  names(m) <- rnam

  deriv <- if(length(object$famfunctions)) eval(object$famfunctions[5]) else dlogis
  dmean.dalpha <- sweep(deriv(lp_full), 2, dvalues, '*')        # n x k
  dmean.dbeta  <- if(p > 0) rowSums(dmean.dalpha) * X else matrix(0., n, 0)
  WD <- dmean.dalpha %*% ytransmat[-1, , drop = FALSE]          # n x myde
  # Claude Sonnet 5 2026-11-20
  # Same flattening order as expandZ()/Zexp elsewhere in this file
  # (column-major by ytrans column): column c = (m-1)*qyde + q, so Z's
  # columns cycle q=1..qyde within each of myde blocks, and WD's
  # columns repeat each m qyde times -- confirmed against expandZ()'s
  # own as.vector(outer(Zrow, ytransrow)) convention.
  dmean.dtau <- Z[, rep(1 : qyde, myde), drop = FALSE] *
                WD[, rep(1 : myde, each = qyde), drop = FALSE]  # n x qm
  dmean.dtheta <- cbind(dmean.dalpha, dmean.dbeta, dmean.dtau)
  mean.var <- diag(dmean.dtheta %*% infoMxop(object$info.matrix, B = t(dmean.dtheta)))
  se <- sqrt(mean.var)
  names(se) <- rnam
  list(mean = m, se.fit = se)
}

##' `Mean()` Closure for a Y-Dependent-Effect `orm()` Fit
##'
##' Internal helper for [Mean.orm()] when `qyde > 0`. [Mean.lrm()]'s own
##' closure takes a bare linear predictor (`lp=`) or design matrix (`X=`)
##' -- sufficient for an ordinary fit, where the mean is a function of a
##' single linear predictor value per observation, but not for a `yde`
##' fit: the Y-dependent effect varies by cutpoint, so computing the mean
##' genuinely needs the `yde` design matrix `Z` as well, at every
##' cutpoint, not just at whichever one a single `lp` was evaluated at --
##' no `lp`-only (or `X`-only) signature can represent that. So, unlike
##' [Mean.lrm()], the function this returns requires `newdata` (the
##' original covariate data frame, from which both `X` and `Z` can be
##' rebuilt) and delegates directly to `predict(object, newdata,
##' type='mean', ...)` -- the answer to "should `Mean.orm` use
##' `predict(type='mean')`": yes, specifically for this case, since
##' `predict()` is the thing that already knows how to build `Z` from
##' `newdata`, which a bare `lp`/`X` cannot supply. Calling the returned
##' function with only `lp=`/`X=` (e.g. via [Predict()]'s own `fun=`
##' mechanism, or [nomogram()]'s direct-lp evaluation, both of which only
##' ever have an `lp`/a single-cutpoint `X` in hand, not `newdata`) stops
##' with an informative error rather than silently using the wrong
##' computation or crashing inside [predict()] on a malformed `newdata`
##' (confirmed directly: passing a numeric design matrix as `newdata=`
##' fails inside `model.frame()`, with a confusing low-level message --
##' this intercepts that case earlier, with a clear one).
##' @param object a `yde`-fitted `orm()` object
##' @param codes see [predict.lrm()]'s own `codes` argument
##' @return a function of `newdata` (required), `conf.int`, and `tmax`
##'   (not yet implemented), paralleling [Mean.lrm()]'s own return value
##'   (a vector, with a `limits` attribute when `conf.int > 0`)
##' @author Frank Harrell
##' @noRd
ormYdeMean <- function(object, codes=FALSE) {
  function(lp=numeric(0), X=numeric(0), newdata=NULL, tmax=NULL, conf.int=0, ...) {
    if(length(tmax))
      stop('tmax is not yet implemented for an orm() fit using yde/ytrans')
    if(! length(newdata))
      stop('Mean() for an orm() fit using yde/ytrans must be called with ',
           'newdata (the original covariate data frame) -- a bare linear ',
           'predictor or design matrix (lp=/X=) cannot represent a ',
           'Y-dependent effect, which varies by level of Y. Call this ',
           'function as m(newdata=...), or use predict(object, newdata, ',
           'type="mean", se.fit=TRUE) directly.')
    se.fit <- conf.int > 0
    r <- predict(object, newdata, type='mean', codes=codes, se.fit=se.fit)
    if(! se.fit) return(r)
    m  <- r$mean
    w  <- qnorm((1 + conf.int) / 2) * r$se.fit
    attr(m, 'limits') <- list(lower = m - w, upper = m + w)
    m
  }
}

##' `ExProb()` Closure for a Y-Dependent-Effect `orm()` Fit
##'
##' Internal helper for [ExProb.orm()] when `qyde > 0`, paralleling
##' [ormYdeMean()]: [ExProb.orm()]'s own closure takes a bare `lp`/`X`,
##' sufficient for an ordinary fit (it gets every cutpoint's probability
##' from one linear predictor via the "shift trick" --
##' `cumprob(sapply(c(1e30, intercepts), '+', lp))` -- valid only because
##' the non-intercept part of `lp` doesn't vary by cutpoint). That shift
##' trick is exactly what breaks for a `yde` fit, so this requires
##' `newdata` and sources the per-cutpoint probabilities from
##' `predict(object, newdata, type='fitted')` instead -- already
##' correctly `yde`-aware (Sec. 14.2 of `orm-yde-implementation.md`) and
##' returning precisely the same `n x k` matrix of `P(Y >= values[j+1])`
##' the shift trick was approximating. `cbind(1, P)` prepends the
##' (always-true) `P(Y >= values[1]) = 1` column, matching the original
##' closure's own `c(1e30, intercepts)` trick for the same purpose. The
##' rest -- the `y=NULL` vs. `y=`-given branching, the step-function
##' `approx(..., method='constant')` lookup for arbitrary queried `y`
##' values, the `'ExProb'`-class return shape -- is copied verbatim from
##' [ExProb.orm()], since none of it is `yde`-specific.
##' @param object a `yde`-fitted `orm()` object
##' @param codes see [ExProb.orm()]'s own `codes` argument
##' @author Frank Harrell
##' @noRd
ormYdeExProb <- function(object, codes=FALSE) {
  function(lp=numeric(0), X=numeric(0), y=NULL, newdata=NULL, conf.int=0, ...) {
    if(! length(newdata))
      stop('ExProb() for an orm() fit using yde/ytrans must be called with ',
           'newdata (the original covariate data frame) -- a bare linear ',
           'predictor or design matrix (lp=/X=) cannot represent a ',
           'Y-dependent effect, which varies by level of Y. Call this ',
           'function as f(newdata=...), or use predict(object, newdata, ',
           'type="fitted") directly.')
    if(conf.int)
      stop('conf.int is not yet implemented for ExProb() on an orm() fit using yde/ytrans')
    values <- if(codes) 1 : length(object$yunique) else as.numeric(object$yunique)
    if(any(is.na(values)))
      stop('values of response levels must be numeric for codes=FALSE')
    yname <- all.vars(object$terms)[1]
    P <- predict(object, newdata, type='fitted')
    if(! is.matrix(P)) P <- matrix(P, nrow = 1)
    prob <- cbind(1, P)
    if(! length(y)) {
      colnames(prob) <- paste('Prob(Y>=', values, ')', sep='')
      result <- structure(list(y=values, prob=prob, yname=yname), class='ExProb')
    } else {
      p <- t(apply(prob, 1,
                 function(probs) {
                   pa <- approx(values, probs, xout=y, f=1, method='constant')$y
                   pa[y < min(values)] <- 1.
                   pa[y > max(values)] <- 0.
                   pa
                 } ) )
      if(length(y) > 1) {
        colnames(p) <- paste('Prob(Y>=', y, ')', sep='')
        p <- list(y=y, prob=p)
      }
      result <- structure(drop(p), class='ExProb')
    }
    result
  }
}

##' `Quantile()` Closure for a Y-Dependent-Effect `orm()` Fit
##'
##' Internal helper for [Quantile.orm()] when `qyde > 0`, paralleling
##' [ormYdeMean()]/[ormYdeExProb()]. [Quantile.orm()]'s own closure
##' builds `cp <- cbind(cumprob(lb), 0)` (or, for `method='discrete'`,
##' `m.cdf <- cbind(1 - cumprob(lb), 1)`) from `lb`, the shift-trick
##' matrix of per-cutpoint linear predictors -- invalid for `yde` for the
##' same reason as [ormYdeExProb()]. Requires `newdata` and substitutes
##' `P <- predict(object, newdata, type='fitted')` for `cumprob(lb)`
##' (verified column-for-column identical in meaning: `P`'s `j`-th column
##' is `P(Y >= values[j+1])`, exactly what `cumprob(lb)[,j]` was). Every
##' other line -- the `ws`/`m.yvals` interpolation-grid construction, the
##' `approx()` inversion, the `'discrete'` method's `which(x >= q)`
##' lookup -- is copied verbatim from [Quantile.orm()], unchanged, since
##' none of it is `yde`-specific; only the source of the cumulative-
##' probability matrix differs.
##' @param object a `yde`-fitted `orm()` object
##' @param codes see [Quantile.orm()]'s own `codes` argument
##' @author Frank Harrell
##' @noRd
ormYdeQuantile <- function(object, codes=FALSE) {
  function(q=numeric(0), lp=numeric(0), X=numeric(0), newdata=NULL, conf.int=0,
           method=c('interpolated', 'discrete'), ...) {
    if(! length(newdata))
      stop('Quantile() for an orm() fit using yde/ytrans must be called with ',
           'newdata (the original covariate data frame) -- a bare linear ',
           'predictor or design matrix (lp=/X=) cannot represent a ',
           'Y-dependent effect, which varies by level of Y. Call this ',
           'function as f(newdata=...), or use predict(object, newdata, ',
           'type="fitted") directly.')
    if(conf.int)
      stop('conf.int is not yet implemented for Quantile() on an orm() fit using yde/ytrans')
    method <- match.arg(method)
    values <- if(codes) 1 : length(object$yunique) else as.numeric(object$yunique)
    if(any(is.na(values)))
      stop('values of response levels must be numeric for codes=FALSE')
    ns <- length(values) - 1L
    trans   <- object$famfunctions
    inverse <- eval(trans[2])
    cumprob <- eval(trans[1])
    P <- predict(object, newdata, type='fitted')
    if(! is.matrix(P)) P <- matrix(P, nrow = 1)
    n <- nrow(P)

    if(method == 'interpolated') {
      m.yvals <- matrix(NA, nrow = n, ncol = ns + 2)
      cp <- cbind(P, 0)
      for(j in 1 : n) {
        ws <- c(0, (cp[j, -ns-1] - cp[j, 1]) / (cp[j, ns] - cp[j, 1]), 1)
        m.yvals[j,] <- (1 - ws) * c(values[1], values) + ws * c(values, values[ns + 1])
      }
      z <- sapply(1 : n,
                 function(i) approx(c(1, cp[i,]), m.yvals[i,],
                                    xout = cumprob(inverse(1 - q)), rule = 2)$y)
    }
    if(method == 'discrete') {
      m.cdf <- cbind(1 - P, 1)
      id <- apply(m.cdf, 1, FUN=function(x) {min(which(x >= q))[1]})
      z <- values[id]
    }
    names(z) <- rownames(P)
    z
  }
}

##' `Survival()` Closure for a Y-Dependent-Effect `orm()` Fit
##'
##' Internal helper for [Survival.orm()] when `qyde > 0`, paralleling
##' [ormYdeMean()]/[ormYdeExProb()]/[ormYdeQuantile()]. [Survival.orm()]'s
##' own closure computes `lpsi <- lps + intercepts[ints]` then
##' `surv <- cump(lpsi)` -- the same shift trick invalid for `yde`
##' (`ints` maps each requested time to a cutpoint index via
##' `approx(values, 1:k, xout=times, method='constant')`, orthogonal to
##' `yde` and reused unchanged here). This instead requires `newdata` and
##' looks the per-(row,time) survival probability up directly from
##' `Pfull <- cbind(predict(object, newdata, type='fitted'), 0)` (an
##' `n x k` matrix -- the appended `0` column matches the original's
##' fake `-Inf` intercept, both meaning "beyond the last modeled
##' cutpoint, survival 0"), via ordinary two-column matrix indexing
##' (`Pfull[cbind(Xrow, ints)]`) rather than `cump(lps + intercepts[ints])`.
##' `surv_last` (this fit's analog of the original's `lpsi_last`, the
##' survival at the last *real* cutpoint, used for the right-censored
##' carry-forward region) is likewise `Pfull[, ns]` rather than
##' `cump(lps + intercepts[k-1])`.
##'
##' Every boundary/censoring-tail branch -- `ts < min_uncens` -> 1,
##' `ts >= max_uncens` -> 0 (no censoring tail), or the tail carry-
##' forward/NA-beyond-follow-up logic (censoring tail present) -- is the
##' same *decision logic* as [Survival.orm()]'s own closure, just setting
##' `surv` directly at those positions instead of going through an `LP`
##' sentinel value (`Inf`/`lpsi_last`) and a final `cump(LP)` -- there is
##' no single scalar "LP" to use as a sentinel for a `yde` fit, since the
##' Y-dependent contribution differs by cutpoint, so the direct-lookup
##' values are set in place instead. `conf.int`/`parallel=TRUE` are not
##' yet implemented (stopped on explicitly) -- a materially larger lift
##' than the point estimate (the original's own SE code inverts parts of
##' `info` keyed by intercept/slope position, which would need the same
##' kind of `tau`-aware generalization as [ormYdeFittedMean()]'s
##' `se.fit=TRUE` path, not attempted in this pass).
##' @param object a `yde`-fitted `orm()` object
##' @param precision see [Survival.orm()]'s own `precision` argument
##' @author Frank Harrell
##' @noRd
ormYdeSurvival <- function(object, precision=7) {
  vals <- as.numeric(object$yunique)
  up <- object$yupper
  if(length(up)) {
    vals <- ifelse(is.infinite(vals), up,
                    ifelse(is.infinite(up), vals, (vals + up) / 2.))
    vals <- as.numeric(vals)
  }
  tl <- object$tail
  if(! length(tl) && length(object$rt_cens_beyond)) {
    rcb <- object$rt_cens_beyond
    tl <- list(type='right', index=length(vals),
               value=rcb$newlevel, tied=NA, range=rcb$range)
  }
  # Claude Sonnet 5 2026-11-20
  # Must NOT write precision=precision as this inner function's own
  # formal default -- that is a self-referencing promise (the formal
  # and the enclosing argument share a name), confirmed directly: it
  # errors "promise already under evaluation" on the first call.
  # Assigning the formal's default via formals()<- afterward, exactly
  # as Survival.orm()'s own non-yde closure does for several of its
  # arguments, evaluates precision in the enclosing scope instead.
  g <- function(times=numeric(0), newdata=NULL, conf.int=0, parallel=FALSE,
                forcedf=FALSE, zero=FALSE, precision=7, ...) {
    if(! length(newdata))
      stop('Survival() for an orm() fit using yde/ytrans must be called with ',
           'newdata (the original covariate data frame) -- a bare linear ',
           'predictor or design matrix (lp=/X=) cannot represent a ',
           'Y-dependent effect, which varies by level of Y. Call this ',
           'function as f(newdata=...), or use predict(object, newdata, ',
           'type="fitted") directly.')
    if(conf.int)
      stop('conf.int is not yet implemented for Survival() on an orm() fit using yde/ytrans')
    if(parallel)
      stop('parallel=TRUE is not yet implemented for Survival() on an orm() fit using yde/ytrans')

    values   <- vals
    tailinfo <- tl
    ranges   <- object$ranges
    len_t    <- length(times)

    if(len_t) {
      mul   <- 10 ^ precision
      times <- round(times * mul) / mul
    }

    zero <- zero && ! len_t && (values[1] != 0)
    if(len_t) times <- sort(unique(times))
    else {
      times <- values
      if(length(tailinfo)) times <- values[- tailinfo$index]
    }

    k <- length(values)

    P <- predict(object, newdata, type='fitted')      # n x (k - 1)
    if(! is.matrix(P)) P <- matrix(P, nrow = 1)
    n  <- nrow(P)
    ns <- ncol(P)
    if(ns != k - 1L)
      stop('logic error in Survival() for an orm() fit using yde/ytrans: ',
           'ncol(fitted)=', ns, ' length(values)-1=', k - 1L)
    Pfull <- cbind(P, 0)                               # n x k; last column: the fake "-Inf" cutpoint

    xrow <- seq_len(n)
    w    <- expand.grid(time = times, Xrow = xrow)
    ts   <- w$time
    ints <- approx(values, 1 : k, xout = ts, method = 'constant', rule = 2)$y

    surv      <- Pfull[cbind(w$Xrow, ints)]
    surv_last <- Pfull[, ns]                            # survival at the last *real* cutpoint, by row

    min_uncens <- ranges$u[1]
    max_uncens <- ranges$u[2]
    surv[ts < min_uncens] <- 1.0
    if(length(tailinfo)) {
      rng    <- tailinfo$range
      max_fu <- rng[2]
      surv[ts > max_fu] <- NA
      use_last <- ts > max_uncens & ts <= max_fu
      surv[use_last] <- surv_last[w$Xrow[use_last]]
    } else surv[ts >= max_uncens] <- 0.0

    w <- data.frame(time = ts, surv, Xrow = w$Xrow)
    if(n == 1) w$Xrow <- NULL
    row.names(w) <- NULL
    if(zero) {
      z <- subset(w, time == min(w$time))
      z$time <- 0; z$surv <- 1
      w <- rbind(z, w)
      row.names(w) <- NULL
    }
    if(forcedf) w
    else if(len_t == 1) w$surv
    else if(n == 1) structure(w$surv, names = paste0('P(T>', w$time, ')'))
    else w
  }
  formals(g)$precision <- precision
  g
}

Mean.lrm <- function(object, codes=FALSE, ...)
{
  ns <- object$non.slopes
  if(ns < 2)
    stop('using this function only makes sense for >2 ordered response categories')
  if(codes) vals <- 1:length(object$freq)
  else {
    vals <- object$yunique
    if(!length(vals)) vals <- names(object$freq)
    vals <- as.numeric(vals)
    if(any(is.na(vals)))
      stop('values of response levels must be numeric for codes=FALSE')
  }
  f <- function(lp=numeric(0), X=numeric(0), tmax=NULL,
                intercepts=numeric(0), slopes=numeric(0),
                info=numeric(0), values=numeric(0),
                interceptRef=integer(0), Ncens, famfunctions=NULL, conf.int=0)
  {
    ns <- length(intercepts)
    lp <- if(length(lp)) lp - intercepts[interceptRef] else matxv(X, slopes) 
    xb <- sapply(intercepts, '+', lp)
    cumprob <- eval(famfunctions[1])
    deriv   <- eval(famfunctions[5])
    P  <- matrix(cumprob(xb), ncol = ns)

    if(! length(tmax)) {
      if(length(Ncens) && sum(Ncens) > 0 && min(1 - P) > 1e-3)
        warning('Computing the mean when the lowest P(Y < y) is ',
                format(min(1 - P)), '\nand tmax omitted will result in only a lower limit to the mean')
    } else {
      if(tmax > max(values))
        stop('tmax=', tmax, ' > maximum observed Y=', format(max(values)))
      values[values > tmax] <- tmax
    }

    P  <- cbind(1, P) - cbind(P, 0)
    m  <- drop(P %*% values)
    names(m) <- names(lp)
    if(conf.int) {
      if(! length(X)) stop('must specify X if conf.int > 0')
      lb <- matrix(sapply(intercepts, '+', lp), ncol = ns)
      dmean.dalpha <- t(apply(deriv(lb),
                              1, FUN=function(x)
                                x * (values[2 : length(values)] - values[1 : ns])))
      dmean.dbeta  <- apply(dmean.dalpha, 1, sum) * X
      dmean.dtheta <- cbind(dmean.dalpha, dmean.dbeta)
      if(getOption('rmsdebug', FALSE)) {prn(infoMxop(info, np=TRUE)); prn(dim(dmean.dtheta))}
      mean.var <- diag(dmean.dtheta %*% infoMxop(info, B=t(dmean.dtheta)))
      w <- qnorm((1 + conf.int) / 2) * sqrt(mean.var)   
      attr(m, 'limits') <- list(lower = m - w, 
                                upper = m + w)
    }
    m
  }
  ## If lrm fit, add information that orm fits have
  family <- object$family
  famf   <- object$famfunctions
  if(! length(family)) {
    family <- 'logistic'
    famf   <- probabilityFamilies$logistic
    }
  ir <- object$interceptRef
  if(!length(ir)) ir <- 1
  # In the following ns is the original number of intercepts
  formals(f) <- list(lp=numeric(0), X=numeric(0), tmax=NULL,
                     intercepts=object$coef[1 : ns],
                     slopes=object$coef[- (1 : ns)],
                     info=object$info.matrix,   
                     values=vals, interceptRef=ir, Ncens=object$Ncens1,
                     famfunctions=famf,
                     conf.int=0)
  f 
}
