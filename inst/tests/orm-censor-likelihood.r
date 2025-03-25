# Check that different configurations of censoring are handled correctly by ormll.f90
# used by orm.fit
# Turn on debugging so that ormll will print subscripts of first terms (ia) and of
# second per-observation terms (ia2) plus the sign used for the first term
# sgn = -1  for Pr(Y = 0) or for left censoring

require(rms)

# Function to access Fortran subroutine used by orm.fit and return
# gradient vector and hessian matrix

rfort <- function(theta, k, y, y2, x=numeric(0), intcens=FALSE, what=3L) {
  m   <- length(theta)
  n   <- length(y)
  p   <- as.integer(m - k)
  nai <- as.integer(if(intcens) 1000000 else 0)
  offset <- numeric(m)
  wt     <- rep(1e0, n)
  penmat <- matrix(0e0, p, p)
  link   <- 1L
  nu     <- 0L
  w <- .Fortran('ormll', n, k, p, x, y, y2, offset, wt, penmat,
                link=link, theta[1:k], theta[-(1:k)],
                logL=numeric(1), grad=numeric(k + p), lpe=numeric(n),
                a=matrix(0e0, (1 - intcens) * k, 2), b=matrix(0e0, p, p), ab=matrix(0e0, k, p),
                intcens, row=integer(nai), col=integer(nai), ai=numeric(nai),
                nai=nai, ne=integer(1),
                urow=integer(nu), ucol=integer(nu), um=numeric(nu), nu=nu, nuu=integer(1),
                what=what, debug=0L, 1L, salloc=integer(1), PACKAGE='rms')
  if(intcens && what == 3L) {
    ne    <- w$ne
    w$a   <- list(row = w$row[1 : ne], col = w$col[1 : ne], a = w$ai[1 : ne])
    w$row <- w$col <- w$ai <- w$nai <- w$ne <- NULL
  }
  info <- infoMxop(w[c('a', 'b', 'ab')])
  z    <- qr(info)
  red  <- z$pivot[-(1 : z$rank)]
  diagzero <- Matrix::diag(info) == 0e0
  if(! identical(which(diagzero), red))
    cat('Redundant rows without zero diagonal\n')
  list(grad=w$grad, info=info, redundant=red)
}

ell <- function(y, y2, long=FALSE) {
  o <- Ocens2ord(Ocens(y, y2))
  k <- length(attr(o, 'levels')) - 1
  if(long) {
    cat('Ocens2ord levels:\n')
    print(attr(o, 'npsurv')$time)
  }

  # Find distinct points whether censored or not
  # u <- unique(sort(c(y[is.finite(y)], y2[is.finite(y2)])))
  u <- unique(sort(c(y, y2)))
  prn(u)
  k2 <- length(u) - 1L
  # Map y and y2 to 0 : k2 preserving censoring
  a <- match(y, u) - 1
  a[is.na(a)] <- -1
  b <- match(y2, u) - 1
  b[is.na(b)] <- k2 + 1
  storage.mode(a) <- 'integer'
  storage.mode(b) <- 'integer'

  init <- qlogis((k2 : 1) / (k2 + 1))
  r <- rfort(init, k2, a, b)
  print(r)
  red <- a %in% r$redundant | b %in% r$redundant
  a[a == -1] <- -Inf
  b[b == k2 + 1] <- Inf
  prn(k2);prn(r$redundant)
  cat('k=', k, '', k2 - length(r$redundant), '\n')

  w <- data.frame(y, y2, A=o[, 1] - 1, B=o[, 2] - 1, a, b,
                  redundant=ifelse(red, '*', ''))
  w
}

# Censoring beyond uncensored values
ell(c(-Inf, 1, 2, 3, 4), c(0, 1, 2, 3, Inf))

# Censoring at outer uncensored values
ell(c(-Inf, 1, 2, 3, 3), c(1, 1, 2, 3, Inf))

# Add interior L and R censored values untied with uncensored values
ell(c(-Inf, 1, -Inf,  2, 2.5, 3, 3), c(1, 1, 1.5, 2, Inf, 3, Inf))

y <- 1:10; y2 <- y
cens <- c(2, 5, 7)
y[cens] <- y[cens] + .01
y2[cens] <- Inf
cbind(y, y2)
ell(y, y2)

y <- 1:10; y2 <- y
cens <- c(2, 3, 5, 7, 9)
y[cens] <- y[cens] + .01
y2[cens] <- Inf
cbind(y, y2)
ell(y, y2)

y <- 1:10; y2 <- y; cens <- 2:4; y2[cens] <- Inf
cbind(y, y2)
ell(y, y2)

y <- 1:10; y2 <- y; cens <- 2:4; y2[cens] <- Inf
y <- c(y, 8); y2 <- c(y2, Inf)
cbind(y, y2)
ell(y, y2)

