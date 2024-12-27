#' Operate on Information Matrices
#'
#' Processes three types of information matrices: ones produced by the `SparseM` package for the `orm` function in `rms` version 6.9-0 and earlier, by the `Matrix` package for version 7.0-0 of `rms`, or plain matrices.  For `Matrix`, the input information matrix is a list with three elements: `a` containing in two columns the diagonal and superdiagonal for intercepts, `b`, a square matrix for the covariates, and `ab` for intercepts x covariates.  If nothing else is specified, the assembled information matrix is returned for `Matrix`, or the original `info` otherwise.  If `p=TRUE`, the number of parameters in the model (number of rows and columns in the whole information matrix) is returned.  If `i` is given, the `i` elements of the inverse of `info` are returned, using efficient calculation to avoid inverting the whole matrix.  Otherwise if `invert=TRUE` or `B` is given without `i`, the efficiently (if `Matrix` or `SparseM`) inverted matrix is returned, or the matrix multiplication of the inverse and `B`.  If both `i` and `B` are given, what is returned is the `i` portion of the inverse of the information matrix, matrix multiplied by `B`.  This is done inside `solve()`.
#' 
#' When inverting `info`, if `info` has a `'scale'` attribute with elements `mean` and `sd`, the scaling is reversed after inverting `info`.
#'
#' @param info an information matrix object
#' @param i integer vector specifying elements returned from the inverse
#' @param invert set to `TRUE` to invert `info` (implied when `i` or `B` is given)
#' @param B multiplier matrix
#' @param np set to `TRUE` to just fetch the total number of parameters (intercepts + betas)
#' @param tol tolerance for matrix inversion singularity
#' @param abort set to `FALSE` to run the `solve` calculation through `try()` without aborting; the user will detect that the operation did not success by examinine `inherits(result, 'try-error')` for being `TRUE`.
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
#' }
infoMxop <- function(info, i, invert=! missing(i) || ! missing(B),
                     B=NULL, np=FALSE, tol=1e-14, abort=TRUE) {
  if(! missing(i) && ! invert)
    stop('i is irrelevant if invert=FALSE')
  
  sc <- attr(info, 'scale')
  
  type <- 'plain'
  if(inherits(info, 'matrix.csr')) type <- 'SparseM'
  else if(! is.matrix(info)) {
    type <- 'Matrix'
    # Assumed to be a 3-element list produced by orm.fit etc.
    a  <- info$a   # intercepts
    b  <- info$b   # betas
    ab <- info$ab  # intercepts x betas
    if(np) return(sum(dim(ab)))
    k  <- nrow(ab) # no. of intercepts = nrow(a)
    p  <- ncol(ab) # no. of betas
    # Simplify if only one intercept, no need for sparseness
    a <- if(k == 1) a[1, 1]
    else Matrix::bandSparse(k, k=c(0,1), diagonals=a, symmetric=TRUE)
    info <- rbind(cbind(a, ab), cbind(t(ab), b))
  }
  if(np) return(ncol(info))

  if(! invert) return(info)
  nv    <- ncol(info)

  # ChatGPT confirmed that extracting submatrices of t(trans) x V x trans equals
  # operating on a submatrix of trans: https://chatgpt.com/share/676e6cb9-bde0-800a-b5f6-0b2c53393ae1
  if(length(sc)) {
    p <- length(sc$mean)
    k <- nv - p
    # t(trans) %*% covariance matrix %*% trans = rescaled cov matrix
    trans <- rbind(cbind(diag(k), matrix(0, nrow=k, ncol=p)),
                   cbind(-matrix(rep(sc$mean / sc$sd, k), ncol=k),
                  diag(1 / as.vector(sc$sd))))
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
  if(missing(i)) {
    v <- if(length(B)) tryit(asm(solv(info, B, tol=tol))) else tryit(asm(solv(info, tol=tol)))
    if(length(sc)) v <- t(trans) %*% v %*% trans
  }
  else {
    # Construct w = a p x r matrix where r = no. desired inverse elements
    # jth column of w has a 1 in i(j) row
    l <- length(i)
    w <- matrix(0., nv, l)
    w[cbind(i, 1 : l)] <- 1
    if(! missing(B)) w <- w %*% B
    v <- tryit(asm(solv(info, w, tol=tol)[i, , drop=FALSE]))
    if(length(sc)) {
      w <- trans[i, i, drop=FALSE]
      v <- t(w) %*% v %*% w
    }
  }
  v
}
