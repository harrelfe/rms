#' Adaptive orm Fit For a Single Continuous Predictor
#'
#' Finds the best fitting `orm` model for a single continuous predictor `x`, starting by finding the link function with the smallest deviance when the predictor is modeled with 4 knots in a restricted cubic spline function.  Then the function finds the number of knots minimizing AIC when this link function is used.  Candidate number of knots are 0 (linear fit), 3, 4, 5, `maxk`.
#'
#' @param x a numeric vector
#' @param y a numeric or factor variable representing an ordinal dependent variable, or an `Ocens` object
#' @param maxk maximum number of knots to try
#' @param ... arguments to orm other than `family`, `x`, `y`
#' @returns the best `orm` fit object
#' @export
#' @md
#' @author Frank Harrell
#'
#' @examples
#' \dontrun{
#' f <- adapt_orm(age, blood_pressure)
#' f$stats['d.f.']   # print no. of parameters for age
#' f$family          # print optimum link found
#' }
adapt_orm <- function(x, y, maxk=6, ...) {
  envi <- .GlobalEnv
  ki   <- min(4, maxk)
  assign('.ki.', ki, envir=envi)
  fi   <- if(ki == 0) orm(y ~ x,            x=TRUE, y=TRUE, ...)
  else                orm(y ~ rcs(x, .ki.), x=TRUE, y=TRUE, ...)
  w    <- Olinks(fi)
  best <- which.min(w[, 'deviance'])
  link <- w[best, 'link']
  if(link != 'logistic')
    fi <- if(ki == 0) orm(y ~ x,            family=link, ...)
      else            orm(y ~ rcs(x, .ki.), family=link, ...)
  fits <- list()
  # Use the best link in finding best # knots
  i  <- 0
  ks <- if(maxk < 3) 0 else c(0, 3 : maxk)
  for(k in ks) {
    assign('.k.', k, envir=envi)   # trick to make modelData work
    f <- if(k == ki) fi
    else if(k == 0) orm(y ~ x,           family=link, ...)
    else            orm(y ~ rcs(x, .k.), family=link, ...)
    if(! f$fail) {
      i <- i + 1
      fits[[i]] <- f
    }
  }
  aic <- sapply(fits, AIC)
  fits[[which.min(aic)]]
}
