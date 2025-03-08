# Check that different configurations of censoring are handled correctly by ormll.f90
# used by orm.fit
# Turn on debugging so that ormll will print subscripts of first terms (ia) and of
# second per-observation terms (ia2) plus the sign used for the first term
# sgn = -1  for Pr(Y = 0) or for left censoring

require(rms)
options(orm.fit.debug=TRUE)

unp <- function(x) {
  i <- grep(paste0('^', x, '$'), ormout)[1] + 1
  x <- ormout[i]
  x <- gsub('  ', ' ', x)
  x <- sub('[1] ', '', x, fixed=TRUE)
  as.numeric(strsplit(x, split=' ')[[1]])
}

w <- function(a, b) {
  y <- Ocens(a, b)
  Y <- Ocens2ord(y, verbose=TRUE)
  f <- attr(Y, 'npsurv')
  attributes(Y) <- attributes(Y)['dim']
  u <- orm.fit(y = y, onlydata=TRUE)
  with(u, cat('\nk=', u$k, ' Ncens=', u$Ncens, '\n'))
  with(f, prn(cbind(time, surv), 'npsurv'))
  ormout <<- capture.output(fit <<- orm.fit(y = y))
  p <- c(1, plogis(coef(fit)))
  prn(cbind(time=fit$yunique, upper=fit$yupper, surv=p), 'orm.fit')
  ia  <- unp('ia')
  ia2 <- unp('ia2')
  sgn <- unp('sgn')
  d <- data.frame(a, b, Y, u$Y)   # ? a, b, y ??
  cat('\n')
  print(d)
  cat('\n')
  print(data.frame(ia, ia2, sgn))
  # cat(z[c(ia, ia + 1, ia2, ia2 + 1, sgn, sgn + 1)], sep='\n')
}

# Censoring beyond uncensored values
w(c(-Inf, 1, 2, 3, 4), c(0, 1, 2, 3, Inf))

# Censoring at outer uncensored values
w(c(-Inf, 1, 2, 3, 3), c(1, 1, 2, 3, Inf))

# Add interior L and R censored values untied with uncensored values
w(c(-Inf, 1, -Inf,  2, 2.5, 3, 3), c(1, 1, 1.5, 2, Inf, 3, Inf))

