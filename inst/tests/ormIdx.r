# Test the exported ormIdx() wrapper around the ormidx Fortran subroutine.
# Compares against an independent pure-R statement of the indexing rules
# documented in ormidx.f90, over every censoring configuration, and also
# checks the edge case of no observation having a second intercept (nb = 0).
require(rms)

ref <- function(y, y2, k) {
  n   <- length(y)
  ia  <- ia2 <- integer(n)
  sgn <- rep(1, n)
  for(i in seq_len(n)) {
    a <- y[i]; b <- y2[i]
    if(a == b) {
      if(a == 0)      { ia[i] <- 1L; sgn[i] <- -1 }
      else if(a == k) ia[i] <- k
      else            { ia[i] <- a; ia2[i] <- a + 1L }
    } else if(a == -1) { ia[i] <- max(b, 1L); sgn[i] <- -1 }
    else if(b > k)                { ia[i] <- min(a + 1L, k) }
    else if(a == 0 && b < k)      { ia[i] <- b + 1L; sgn[i] <- -1 }
    else if(a > 0 && b == k)      { ia[i] <- a }
    else if(a > 0 && b < k && a < b) { ia[i] <- a; ia2[i] <- b + 1L }
    else stop('unhandled')
  }
  list(ia=ia, ia2=ia2, sgn=sgn, ib=which(ia2 > 0), nb=sum(ia2 > 0))
}

k <- 6L
# uncensored 0..k, left censored (-1, b), right censored (a, k+1),
# interval censored a < b
# (an interval spanning all categories, a=0 and b=k, is uninformative and is
# not a configuration ormidx handles, so it is not included)
y  <- c(0:k,  -1L, -1L, -1L,  0L,  3L,  k,  0L,  2L, 1L,  2L)
y2 <- c(0:k,   0L,  3L,  k,   k + 1L, k + 1L, k + 1L, 4L, k, 2L, 5L)

w <- ormIdx(y, y2, k)
r <- ref(y, y2, k)
for(nm in c('ia', 'ia2', 'sgn', 'ib', 'nb'))
  stopifnot(isTRUE(all.equal(as.numeric(w[[nm]]), as.numeric(r[[nm]]))))
stopifnot(is.integer(w$ia), is.integer(w$ia2), is.integer(w$ib),
          is.numeric(w$sgn), all(w$ia >= 1L & w$ia <= k),
          all(w$ia2 >= 0L & w$ia2 <= k))
# ib must be exactly the observations with a second intercept
stopifnot(identical(w$ib, which(w$ia2 > 0L)), w$nb == length(w$ib))

# default y2 = y (no censoring)
stopifnot(identical(ormIdx(0:k, k=k)$ia, ormIdx(0:k, 0:k, k)$ia))

# nb = 0: seq_len() must give an empty ib, not a bogus element
z <- ormIdx(c(0L, k), c(0L, k), k)
stopifnot(z$nb == 0L, length(z$ib) == 0L)

# error handling
stopifnot(inherits(try(ormIdx(3L, 1L, k), silent=TRUE), 'try-error'))   # y > y2
stopifnot(inherits(try(ormIdx(c(1L, NA), c(1L, 1L), k), silent=TRUE), 'try-error'))
stopifnot(inherits(try(ormIdx(1.5, 1.5, k), silent=TRUE), 'try-error'))
cat('ormIdx tests passed\n')
