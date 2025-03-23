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
  list(grad=w$grad, info=info, redundant=red)
}

ell <- function(y, y2) {
  # Find distinct points whether censored or not
  u <- unique(sort(c(y[is.finite(y)], y2[is.finite(y2)])))
  k <- length(u) - 1
  # Map y and y2 to 0 : k preserving censoring
  a <- match(y, u) - 1
  a[is.na(a)] <- -1
  b <- match(y2, u) - 1
  b[is.na(b)] <- k + 1
  storage.mode(a) <- 'integer'
  storage.mode(b) <- 'integer'
  cat('k=', k, '\n')

  init <- qlogis((k : 1) / (k + 1))
  r <- rfort(init, k, a, b)
  red <- a %in% r$redundant | b %in% r$redundant
  print(data.frame(a, b, redundant=ifelse(red, '*', '')))
  r
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

