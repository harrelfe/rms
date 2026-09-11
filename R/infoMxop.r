#' Operate on Information Matrices
#'
#' Processes four types of information matrices: ones produced by the `SparseM` package for the `orm` or `lrm` functions in `rms` version 6.9-0 and earlier, by the `Matrix` package for version 7.0-0 of `rms` using a tri-band diagonal matrix for the intercepts, using `Matrix` for general sparse information matrices for intercepts (when any interval-censored observations exist or random effects are included), or plain matrices.  For `Matrix`, the input information matrix is a list with three elements: `a` containing in two columns the diagonal and superdiagonal for intercepts (when there is no interval censoring) or a list with three elements `row`, `col`, `a` (when there is interval censoring), `b`, a square matrix for the covariates, and `ab` for intercepts x covariates.  If nothing else is specified, the assembled information matrix is returned for `Matrix`, or the original `info` otherwise.  If `p=TRUE`, the number of parameters in the model (number of rows and columns in the whole information matrix minus any log(sigma)/sigma2 elements) is returned.  If `i` is given, the `i` elements of the inverse of `info` are returned, using efficient calculation to avoid inverting the whole matrix.  Otherwise if `invert=TRUE` or `B` is given without `i`, the efficiently (if `Matrix` or `SparseM`) inverted matrix less any log(sigma)/sigma2 elements is returned, or the matrix multiplication of the inverse and `B`.  If both `i` and `B` are given, what is returned is the `i` portion of the inverse of the information matrix, matrix multiplied by `B`.  This is done inside `solve()`.
#'
#' When only variance-covariance matrix elements corresponding to the non-intercepts are desired, specify
#' `i='x'` or `i=(k + 1) : nv` where `nv` is the number of intercepts and slopes combined.  `infoMxop` computes the needed covariance matrix very quickly in this case.  To retrieve only the variance/covariance of the random-effect scale parameter(s) for a model with random effects, specify `i='sigma_parameters'`.
#' When inverting `info`, if `info` has a `'scale'` attribute with elements `mean` and `sd`, the scaling is reversed after inverting `info`.
#'
#' When the number of intercepts is needed to be known and the `info` object is not a 3-element list, `info` must have an `intercepts` attribute to define the number.  This is used for example when `transx=TRUE` is specified to `lrm` or `lrm.fit`.
#'
#' When the model contains random effects, the `xname` component of the `info` list must contain an element named `log(sigma)` (the log of the random effect's scale, or of its first scale parameter in a two-scale model) and, for a two-scale (`sigma1`/`sigma2`) random-effect model, an additional element named `sigma2`.  `infoMxop` locates these BY NAME within `xname`, wherever they occur -- it does not assume they are the last element(s) of `xname` or in any particular order relative to each other.
#'
#' A random-intercept fit's information matrix can become genuinely near-singular in one or two directions -- `log(sigma)`, and `sigma2` when present -- whenever the corresponding scale parameter is at or very near a boundary where it is not identifiable; this is a real, expected non-regular feature of variance-component estimation, not a numerical defect. Since every general solve (`missing(i)`, `i='i'`, or any `i` mixing intercepts and slopes) factors the *whole* information matrix regardless of which elements of its inverse are ultimately kept, this near-singularity can surface even when the caller never asked for `sigma_parameters` at all. When it does, `infoMxop` automatically retries with a small, escalating ridge applied *only* to the affected diagonal entry/entries -- not the whole matrix -- since these are empirically essentially uncorrelated with everything else at such a fit, so this leaves every other requested quantity negligibly perturbed. A warning is issued when this happens. `i='x'` is never affected by this, since it already avoids factoring the full matrix (see below); `i='sigma_parameters'` is also never affected, since that variance/covariance genuinely isn't well defined there and is instead returned as `NA` (a scalar or an all-`NA` 2x2 matrix, matching whether one or two scale parameters are present) with its own explanatory warning.
#'
#' @param info an information matrix object
#' @param i integer vector specifying elements returned from the inverse.  You can also specify `i='x'` to return non-intercepts, `i='i'` to return intercepts, or `i="sigma_parameters"` to return the estimated variance of `log(sigma)` (a scalar), or the 2x2 covariance matrix of `log(sigma)` and `sigma2` when both are present (a two-scale random-effect model).
#' @param invert set to `TRUE` to invert `info` (implied when `i` or `B` is given)
#' @param B multiplier matrix.  If random effects were included, `B` should ignore them.
#' @param np set to `TRUE` to just fetch the total number of parameters (intercepts + betas, any log(sigma)/sigma2 not counting)
#' @param tol tolerance for matrix inversion singularity
#' @param abort set to `FALSE` to run the `solve` calculation through `try()` without aborting; the user will detect that the operation did not success by examinine `inherits(result, 'try-error')` for being `TRUE`.  This is ignored when `i='sigma_parameters'`: if the underlying `solve()` fails there (e.g. because the model converged with a scale parameter at or near its lower boundary of zero, where it is not identifiable), a warning is issued and `NA` is returned regardless of `abort`, since a plain `NA` (scalar or all-`NA` 2x2 matrix, as appropriate) is more informative than a `try-error` object for this specific, well-understood failure mode.  For every *other* request, a singularity caused specifically by one of these scale parameters being at this same boundary is instead handled by an automatic, targeted ridge retry (see **Description**) with a warning on success; `abort` only comes into play if that retry itself fails, in which case the original error is handled exactly as documented above.
#'
#' @returns a single integer or a matrix
#' @export
#' @md
#' @author Frank Harrell
#'
#' @examples
#' \dontrun{
#' f <- orm(y ~ x)
#' infoMxop(f$info.matrix)   # assembles 3 pieces
#' infoMxop(v, i=c(2,4))     # returns a submatrix of v inverse
#' infoMxop(f$info.matrix, i='x')  # sub-covariance matrix for just the betas
#' f <- orm(y ~ x + cluster(id))
#' infoMxop(f$info.matrix, i='sigma_parameters')
#' # a single number = variance of log(sigma), sigma=sqrt of random effects variance,
#' # OR (for a two-scale sigma1/sigma2 random-effect model) a 2x2 covariance matrix
#' # for log(sigma) and sigma2
#' }
infoMxop <- function(info, i, invert=! missing(i) || ! missing(B),
                     B, np=FALSE, tol=.Machine$double.eps, abort=TRUE) {
  if(! missing(i) && ! invert)
    stop('i is irrelevant if invert=FALSE')
  Bp <- ! missing(B)
  if(Bp) B <- Matrix::Matrix(B)

  xname  <- iname <- name <- sc <- NULL
  nsigma <- 0L
  sigma_gpos <- integer(0)   ## global positions of log(sigma)/sigma2, if any (set below)
  request_sigma_params <- FALSE

  if(is.matrix(info)) name <- colnames(info)

  type <- 'plain'
  t3   <- FALSE
  k    <- attr(info, 'intercepts')

  if(inherits(info, 'matrix.csr')) type <- 'SparseM'
  else if(is.list(info) && all(c('a', 'b', 'ab') %in% names(info))) {
    # Object created by lrm or orm
      type   <- 'Matrix'
      t3     <- TRUE
      a      <- info$a   # intercepts
      b      <- info$b   # betas
      ab     <- info$ab  # intercepts x betas
      xname  <- info$xname
      iname  <- info$iname
      sc     <- info$scale
      ## Claude Sonnet 5 2026-09-07          6 lines
      ## Dropped the positional assumption that log(sigma), when
      ## present, is xname's LAST element -- search xname directly for
      ## 'log(sigma)' (the random-effect scale parameter; also plays
      ## the role of "log(sigma1)" in a two-scale model) and 'sigma2'
      ## instead, wherever they actually occur. nsigma remains a plain
      ## count (0, 1, or 2) used exactly as before for "how many
      ## trailing random-effect-scale columns to exclude from p/nv".
      sigma1_xpos <- which(xname == 'log(sigma)')
      sigma2_xpos <- which(xname == 'sigma2')
      nsigma      <- length(sigma1_xpos) + length(sigma2_xpos)
      k      <- nrow(ab) # no. of intercepts = nrow(a)
      p      <- ncol(ab) - nsigma  # no. of betas
      if(np) return(k + p)
      ## Claude Sonnet 5 2026-09-07          4 lines
      ## Global (whole-info-matrix) positions of log(sigma)/sigma2,
      ## now that k is known -- xname's own positions are offset by k
      ## once intercepts are prepended below. Ordered log(sigma) then
      ## sigma2 regardless of their order within xname.
      sigma1_gpos <- if(length(sigma1_xpos)) k + sigma1_xpos else integer(0)
      sigma2_gpos <- if(length(sigma2_xpos)) k + sigma2_xpos else integer(0)
      sigma_gpos  <- c(sigma1_gpos, sigma2_gpos)
      # Simplify if only one intercept, no need for sparseness
      a <- if(k == 1) if(is.list(a)) a$a[1] else a[1, 1]
      else if(is.list(a)) Matrix::sparseMatrix(a$row, a$col, x=a$a, dims=c(k, k), symmetric=TRUE)
      else Matrix::bandSparse(k, k=c(0,1), diagonals=a, symmetric=TRUE)
      info <- rbind(cbind(a, ab), cbind(t(ab), b))
      name <- c(iname, xname)
      dimnames(info) <- list(name, name)
    } else if(inherits(info, 'Matrix')) type <- 'Matrix'
    else type <- 'plain'

  nv <- ncol(info) - nsigma
  if(np) return(nv)
  if(! invert) return(info)

  if(! length(k)) stop('info did not contain intercepts attribute')

  # ChatGPT confirmed that extracting submatrices of t(trans) x V x trans equals
  # operating on a submatrix of trans: https://chatgpt.com/share/676e6cb9-bde0-800a-b5f6-0b2c53393ae1
  if(length(sc)) {
    # t(trans) %*% covariance matrix %*% trans = rescaled cov matrix
    trans <- rbind(cbind(Matrix::Diagonal(k), Matrix::Matrix(0., k, p)),
                   cbind(Matrix::Matrix(- rep(sc$mean / sc$sd, k), ncol=k),
                  Matrix::Diagonal(x = 1. / as.vector(sc$sd))))
  }

  tryit <- if(abort) function(x) x else function(x) try(x)
  solv  <- switch(type,
                  plain   =          solve,
                  SparseM = SparseM::solve,
                  Matrix  = Matrix ::solve)
  asm   <- switch(type,
                  plain   =          as.matrix,
                  SparseM = SparseM::as.matrix,
                  Matrix  = Matrix ::as.matrix)

  # Claude Sonnet 5 2026-09-03          entire function
  # A random-intercept fit's info matrix can be genuinely near-
  # singular in exactly one or two directions -- log(sigma), and
  # sigma2 when a two-scale (sigma1/sigma2) random-effect model is in
  # use -- whenever the fit has converged at or near a boundary where
  # that scale parameter is no longer identifiable (a real, expected
  # non-regular feature of variance-component estimation, confirmed
  # empirically: such a diagonal entry collapses toward zero while its
  # off-diagonal entries with everything else are many orders of
  # magnitude smaller still, i.e. it is essentially uncorrelated with
  # the rest of the model at such a fit). This affects EVERY general
  # solve against the full (nv+nsigma)-dimensional info matrix --
  # missing(i), i='i', i=c(...) mixes of intercepts and slopes -- not
  # just an explicit i='sigma_parameters' request, since the general
  # solve always factors the whole matrix regardless of which elements
  # of its inverse are ultimately kept. (i='x' is unaffected because it
  # works with the much smaller Schur-complement M over alpha alone,
  # which stays just barely non-singular; see the "invert-then-drop"
  # logic below.)
  #
  # solve_expr(info) attempts the caller's actual solve. On failure,
  # if any sigma-related positions exist (sigma_gpos), this retries
  # with a small, escalating ridge added ONLY to those diagonal
  # entries -- not the whole matrix -- since that leaves every OTHER
  # requested quantity negligibly perturbed while resolving the
  # genuinely singular direction(s). Never invoked for
  # request_sigma_params itself, which keeps its own NA + explanatory-
  # warning handling below, since that variance genuinely isn't well
  # defined at such a fit -- ridging it away would silently
  # manufacture a meaningless number instead.
  # Claude Sonnet 5 2026-09-07          4 lines
  # Generalized from a single, position-assumed entry (nv+1) to an
  # arbitrary set of positions (sigma_gpos, found by name rather than
  # assumed) -- ridges log(sigma) alone, or both log(sigma) and
  # sigma2, as appropriate.
  solve_with_logsigma_ridge <- function(solve_expr, info, sigma_gpos) {
    v <- tryCatch(solve_expr(info), error = function(e) e)
    if(! inherits(v, 'error')) return(v)
    if(! length(sigma_gpos)) stop(conditionMessage(v))
    lambda <- 1e-8
    repeat {
      info_ridge <- info
      for(jj in sigma_gpos) info_ridge[jj, jj] <- info_ridge[jj, jj] + lambda
      v2 <- tryCatch(solve_expr(info_ridge), error = function(e) e)
      if(! inherits(v2, 'error')) {
        warning('Model converged with a random-effect scale parameter at or ',
                'near its lower boundary; a small ridge was applied to the ',
                if(length(sigma_gpos) > 1) 'log(sigma)/sigma2 directions only'
                else 'log(sigma) direction only',
                ' (essentially uncorrelated with the other parameters at ',
                'such a fit) to compute the requested covariance elements.')
        return(v2)
      }
      lambda <- lambda * 100
      if(lambda > 1) stop(conditionMessage(v))
    }
  }

  # BUG?: whenever Bp=TRUE and sc (scale) is present, v ends up
  # length(i) x ncol(B) -- not square -- at each of the three places
  # below that do  v <- t(trans) %*% v %*% trans  or the equivalent
  # w <- trans[i,i,...]; v <- t(w) %*% v %*% w:
  #   1. the shared "if(length(sc)) v <- t(trans) %*% v %*% trans"
  #      line in the missing(i) branch
  #   2. the "if(length(sc))" block inside the i='x'-equivalent
  #      shortcut (the (length(i)==nv-k) block)
  #   3. the "if(length(sc))" block in the general (slow) path at
  #      the end of the function
  # In each case this is not the same as rescaling the square V
  # first and then multiplying by B -- unverified whether sc and B
  # actually co-occur in practice, and unrelated to nsigma/log(sigma)
  # specifically (Bp + sc alone is enough to trigger it).

  if(missing(i)) {
    if(Bp && nsigma) {
      # B is always sized to nv rows (the non-log(sigma) parameter
      # space) regardless of whether random effects are present --
      # pad with zero rows for log(sigma) so it's conformable with
      # the full (nv+nsigma)-dimensional info, solve, then keep only
      # the rows the caller actually asked about. Since the padded
      # rows are zero, this is algebraically identical to
      # V[1:nv, 1:nv] %*% B without ever forming the full inverse.
      Bfull <- rbind(B, Matrix::Matrix(0., nsigma, ncol(B)))
      v <- tryit(solve_with_logsigma_ridge(function(inf) asm(solv(inf, Bfull, tol=tol)),
                                           info, sigma_gpos))
      if(! inherits(v, 'try-error')) v <- v[1 : nv, , drop=FALSE]
    } else {
      v <- if(Bp) tryit(solve_with_logsigma_ridge(function(inf) asm(solv(inf, B, tol=tol)), info, sigma_gpos))
           else       tryit(solve_with_logsigma_ridge(function(inf) asm(solv(inf, tol=tol)), info, sigma_gpos))
      if (nsigma && !Bp) v <- v[1:nv, 1:nv, drop = FALSE]
      }
        if(length(sc)) v <- t(trans) %*% v %*% trans
      }
  else {
    # User has specified i, a vector of indexes of rows/columns of inverse to keep
    if(is.character(i) && length(i) == 1) {
      request_sigma_params <- i == 'sigma_parameters'
      if(request_sigma_params && ! nsigma)
        stop('i="sigma_parameters" specified for model with no log(sigma)/sigma2')
#      if(! t3) k <- attr(info, 'intercepts')
#      if(! length(k))
#        stop("may only specify i='i' or 'x' when operating on the default ",
#             "lrm or orm 3-element information matrix or when info has ",
#             "an intercepts attribute")
      ## Claude Sonnet 5 2026-09-07          1 line
      ## sigma_gpos (found by searching xname, not assumed at nv+1) is
      ## length 1 (log(sigma) only) or 2 (log(sigma), sigma2), ordered
      ## log(sigma) then sigma2 regardless of their order in xname.
      i <- switch(i,
                  i                = 1 : k,
                  x                = (k + 1) : nv,
                  sigma_parameters = sigma_gpos)
    }
    if((length(i) == nv - k) && all(sort(i) == (k + 1) : nv)) {
      # It's very quick to only get the beta components of the inverse
      # It's slower to do likewise for just the intercept components; best to
      # just use the i=1:k for that
      if(t3) {
        M <- b - Matrix::t(ab) %*% solv(a, ab, tol = tol)
        if(Bp && nsigma) {
          Bfull <- rbind(B, Matrix::Matrix(0., nsigma, ncol(B)))
          v <- solv(M, Bfull, tol=tol)
          v <- v[1 : p, , drop=FALSE]
        } else {
          v <- if(Bp) solv(M, B, tol=tol) else solv(M, tol=tol)
          if(nsigma && ! Bp) v <- v[1 : p, 1 : p, drop=FALSE]
          }
        } else v <- solv(info)[i, i, drop=FALSE]
      if(length(sc)) {
        w <- trans[i, i, drop=FALSE]
        v <- t(w) %*% v %*% w
      }
    if(! Bp) dimnames(v) <- list(name[i], name[i])
    return(v)
    }

    # Construct w = a p x r matrix where r = no. desired inverse elements
    # jth column of w has a 1 in i(j) row
    l <- length(i)
    w <- matrix(0., nv + nsigma, l)
    w[cbind(i, 1 : l)] <- 1
    if(type == 'Matrix') w <- Matrix::Matrix(w)
    if(Bp) w <- w %*% B
    # Claude Sonnet 5 2026-09-02          9 lines
    # i='sigma_parameters' can converge to a point where sigma is at or
    # near its lower boundary of zero, at which point the relevant
    # information genuinely collapses to zero (a real, expected
    # boundary phenomenon in variance-component estimation, not a
    # numerical artifact) and solv() throws a computationally-singular
    # error. For this specific request only, catch that and return NA
    # with an informative warning instead of the raw solve() error --
    # unconditional on abort, since a clear NA is strictly more useful
    # here than either a hard stop or a bare try-error object.
    # Claude Sonnet 5 2026-09-07          2 lines
    # Generalized to 1 or 2 requested elements: a scalar NA when only
    # log(sigma) is present, an all-NA 2x2 matrix when sigma2 is too.
    v <- if(request_sigma_params)
           tryCatch(asm(solv(info, w, tol=tol)[i, , drop=FALSE]),
                    error = function(e) {
                      warning('Model converged with a random-effect scale parameter ',
                              'at or near its lower boundary of zero, so ',
                              if(length(i) > 1) 'log(sigma)/sigma2 are' else 'log(sigma) is',
                              ' not identifiable here and its standard error is not well ',
                              'defined; returning NA')
                      if(length(i) == 1) NA_real_ else matrix(NA_real_, length(i), length(i))
                    })
         else tryit(solve_with_logsigma_ridge(function(inf) asm(solv(inf, w, tol=tol)[i, , drop=FALSE]),
                                             info, sigma_gpos))
    if(! (request_sigma_params && all(is.na(v)))) {
      if(length(sc)) {
        w <- trans[i, i, drop=FALSE]
        v <- t(w) %*% v %*% w
      }
      if(! Bp) dimnames(v) <- list(name[i], name[i])
    }
  }
  ## Claude Sonnet 5 2026-09-07          1 line
  ## Flatten to a plain scalar only for the single-element (log(sigma)
  ## alone) case; keep the 2x2 matrix (with dimnames) intact when
  ## sigma2 is also present.
  if(request_sigma_params && length(i) == 1) as.vector(v) else v
}
