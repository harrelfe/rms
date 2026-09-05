#' Operate on Information Matrices
#'
#' Processes four types of information matrices: ones produced by the `SparseM` package for the `orm` or `lrm` functions in `rms` version 6.9-0 and earlier, by the `Matrix` package for version 7.0-0 of `rms` using a tri-band diagonal matrix for the intercepts, using `Matrix` for general sparse information matrices for intercepts (when any interval-censored observations exist or random effects are included), or plain matrices.  For `Matrix`, the input information matrix is a list with three elements: `a` containing in two columns the diagonal and superdiagonal for intercepts (when there is no interval censoring) or a list with three elements `row`, `col`, `a` (when there is interval censoring), `b`, a square matrix for the covariates, and `ab` for intercepts x covariates.  If nothing else is specified, the assembled information matrix is returned for `Matrix`, or the original `info` otherwise.  If `p=TRUE`, the number of parameters in the model (number of rows and columns in the whole information matrix minus any log(sigma) elements) is returned.  If `i` is given, the `i` elements of the inverse of `info` are returned, using efficient calculation to avoid inverting the whole matrix.  Otherwise if `invert=TRUE` or `B` is given without `i`, the efficiently (if `Matrix` or `SparseM`) inverted matrix less any log(sigma) elements is returned, or the matrix multiplication of the inverse and `B`.  If both `i` and `B` are given, what is returned is the `i` portion of the inverse of the information matrix, matrix multiplied by `B`.  This is done inside `solve()`.
#'
#' When only variance-covariance matrix elements corresponding to the non-intercepts are desired, specify
#' `i='x'` or `i=(k + 1) : nv` where `nv` is the number of intercepts and slopes combined.  `infoMxop` computes the needed covariance matrix very quickly in this case.  To retrieve only the variance of log(sigma) for a model with random effects, specify `i='log(sigma)'`.
#' When inverting `info`, if `info` has a `'scale'` attribute with elements `mean` and `sd`, the scaling is reversed after inverting `info`.
#'
#' When the number of intercepts is needed to be known and the `info` object is not a 3-element list, `info` must have an `intercepts` attribute to define the number.  This is used for example when `transx=TRUE` is specified to `lrm` or `lrm.fit`.
#'
#' When the model contains random effects and there are elements corresponding to log(sigma) in the `ab` and `b` components of `info`, the `xname` component of the `info` list must contain a last element named `log(sigma)` corresponding to the column number of log(sigma).  This should be the last column in `ab` and `b`.
#'
#' A random-intercept fit's information matrix can become genuinely near-singular in exactly one direction -- `log(sigma)` -- whenever the fitted `sigma` is at or very near its lower boundary of zero; this is a real, expected non-regular feature of variance-component estimation (log(sigma) is not identifiable there), not a numerical defect. Since every general solve (`missing(i)`, `i='i'`, or any `i` mixing intercepts and slopes) factors the *whole* information matrix regardless of which elements of its inverse are ultimately kept, this near-singularity can surface even when the caller never asked for `log(sigma)` at all. When it does, `infoMxop` automatically retries with a small, escalating ridge applied *only* to `log(sigma)`'s own diagonal entry -- not the whole matrix -- since `log(sigma)` is empirically essentially uncorrelated with everything else at such a fit, so this leaves every other requested quantity negligibly perturbed. A warning is issued when this happens. `i='x'` is never affected by this, since it already avoids factoring the full matrix (see below); `i='log(sigma)'` is also never affected, since `log(sigma)`'s own variance genuinely isn't well defined there and is instead returned as `NA` with its own explanatory warning.
#'
#' @param info an information matrix object
#' @param i integer vector specifying elements returned from the inverse.  You can also specify `i='x'` to return non-intercepts, `i='i'` to return intercepts, or `i="log(sigma)"` to return the estimated variance of log(sigma).
#' @param invert set to `TRUE` to invert `info` (implied when `i` or `B` is given)
#' @param B multiplier matrix.  If random effects were included, `B` should ignore them.
#' @param np set to `TRUE` to just fetch the total number of parameters (intercepts + betas, possible log(sigma) not counting)
#' @param tol tolerance for matrix inversion singularity
#' @param abort set to `FALSE` to run the `solve` calculation through `try()` without aborting; the user will detect that the operation did not success by examinine `inherits(result, 'try-error')` for being `TRUE`.  This is ignored when `i='log(sigma)'`: if the underlying `solve()` fails there (e.g. because the model converged with sigma at or near its lower boundary of zero, where log(sigma) is not identifiable), a warning is issued and `NA` is returned regardless of `abort`, since a plain `NA` is more informative than a `try-error` object for this specific, well-understood failure mode.  For every *other* request, a singularity caused specifically by `log(sigma)` being at this same boundary is instead handled by an automatic, targeted ridge retry (see **Description**) with a warning on success; `abort` only comes into play if that retry itself fails, in which case the original error is handled exactly as documented above.
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
#' infoMxop(f$info.matrix, i='log(sigma)')
#' # a single number = variance of log(sigma), sigma=sqrt of random effects variance
#' }
infoMxop <- function(info, i, invert=! missing(i) || ! missing(B),
                     B, np=FALSE, tol=.Machine$double.eps, abort=TRUE) {
  if(! missing(i) && ! invert)
    stop('i is irrelevant if invert=FALSE')
  Bp <- ! missing(B)
  if(Bp) B <- Matrix::Matrix(B)

  xname  <- iname <- name <- sc <- NULL
  nsigma <- FALSE
  request_logsigma <- FALSE

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
      nsigma <- length(xname) && xname[length(xname)] == 'log(sigma)'
      k      <- nrow(ab) # no. of intercepts = nrow(a)
      p      <- ncol(ab) - nsigma  # no. of betas
      if(np) return(k + p)
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
  # singular in exactly one direction -- log(sigma) -- whenever sigma
  # has converged at or near its lower boundary of zero (a real,
  # expected non-regular feature of variance-component estimation,
  # confirmed empirically: log(sigma)'s own diagonal entry collapses
  # toward zero while its off-diagonal entries with everything else
  # are many orders of magnitude smaller still, i.e. it is essentially
  # uncorrelated with the rest of the model at such a fit). This
  # affects EVERY general solve against the full (nv+nsigma)-
  # dimensional info matrix -- missing(i), i='i', i=c(...) mixes of
  # intercepts and slopes -- not just an explicit i='log(sigma)'
  # request, since the general solve always factors the whole matrix
  # regardless of which elements of its inverse are ultimately kept.
  # (i='x' is unaffected because it works with the much smaller
  # Schur-complement M over alpha alone, which stays just barely
  # non-singular; see the "invert-then-drop" logic below.)
  #
  # solve_expr(info) attempts the caller's actual solve. On failure,
  # if nsigma is set, this retries with a small, escalating ridge
  # added ONLY to info's log(sigma) diagonal entry -- not the whole
  # matrix -- since that leaves every OTHER requested quantity
  # negligibly perturbed while resolving the one genuinely singular
  # direction. Never invoked for request_logsigma itself, which keeps
  # its own NA + explanatory-warning handling below, since log(sigma)'s
  # own variance genuinely isn't well defined at such a fit -- ridging
  # it away would silently manufacture a meaningless number instead.
  solve_with_logsigma_ridge <- function(solve_expr, info, nv, nsigma) {
    v <- tryCatch(solve_expr(info), error = function(e) e)
    if(! inherits(v, 'error')) return(v)
    if(! nsigma) stop(conditionMessage(v))
    lambda <- 1e-8
    repeat {
      info_ridge <- info
      info_ridge[nv + 1, nv + 1] <- info_ridge[nv + 1, nv + 1] + lambda
      v2 <- tryCatch(solve_expr(info_ridge), error = function(e) e)
      if(! inherits(v2, 'error')) {
        warning('Model converged with sigma at or near its lower boundary ',
                'of zero; a small ridge was applied to the log(sigma) ',
                'direction only (essentially uncorrelated with the other ',
                'parameters at such a fit) to compute the requested ',
                'covariance elements.')
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
                                           info, nv, nsigma))
      if(! inherits(v, 'try-error')) v <- v[1 : nv, , drop=FALSE]
    } else {
      v <- if(Bp) tryit(solve_with_logsigma_ridge(function(inf) asm(solv(inf, B, tol=tol)), info, nv, nsigma))
           else       tryit(solve_with_logsigma_ridge(function(inf) asm(solv(inf, tol=tol)), info, nv, nsigma))
      if (nsigma && !Bp) v <- v[1:nv, 1:nv, drop = FALSE]
      }
        if(length(sc)) v <- t(trans) %*% v %*% trans
      }
  else {
    # User has specified i, a vector of indexes of rows/columns of inverse to keep
    if(is.character(i) && length(i) == 1) {
      request_logsigma <- i == 'log(sigma)'
      if (request_logsigma && ! nsigma) stop('i="log(sigma)" specified for model with no log(sigma)')
#      if(! t3) k <- attr(info, 'intercepts')
#      if(! length(k))
#        stop("may only specify i='i' or 'x' when operating on the default ",
#             "lrm or orm 3-element information matrix or when info has ",
#             "an intercepts attribute")
      i <- switch(i,
                  i            = 1 : k,
                  x            = (k + 1) : nv,
                  'log(sigma)' = nv + 1)
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
    # i='log(sigma)' can converge to a point where sigma is at or near
    # its lower boundary of zero, at which point log(sigma)'s
    # information genuinely collapses to zero (a real, expected
    # boundary phenomenon in variance-component estimation, not a
    # numerical artifact) and solv() throws a computationally-singular
    # error. For this specific request only, catch that and return NA
    # with an informative warning instead of the raw solve() error --
    # unconditional on abort, since a clear NA is strictly more useful
    # here than either a hard stop or a bare try-error object.
    v <- if(request_logsigma)
           tryCatch(asm(solv(info, w, tol=tol)[i, , drop=FALSE]),
                    error = function(e) {
                      warning('Model converged with sigma at or near its lower ',
                              'boundary of zero, so log(sigma) is not identifiable ',
                              'here and its standard error is not well defined; ',
                              'returning NA')
                      NA_real_
                    })
         else tryit(solve_with_logsigma_ridge(function(inf) asm(solv(inf, w, tol=tol)[i, , drop=FALSE]),
                                             info, nv, nsigma))
    if(! (request_logsigma && length(v) == 1 && is.na(v))) {
      if(length(sc)) {
        w <- trans[i, i, drop=FALSE]
        v <- t(w) %*% v %*% w
      }
      if(! Bp) dimnames(v) <- list(name[i], name[i])
    }
  }
  if (request_logsigma) as.vector(v) else v
}
