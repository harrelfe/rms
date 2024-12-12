#' Create Ordinal Variables With a Given Precision
#'
#' For a factor variable `y`, uses existing factor levels and codes the output `y` as integer.  For a character `y`, converts to `factor` and does the same.  For a numeric `y` that is integer, leaves the levels intact and codes `y` as consecutive positive integers corresponding to distinct values in the data.  For numeric `y` that contains any non-integer values, rounds `y` to `precision` decimal places to the right before finding the distinct values.
#'
#' This function is used to prepare ordinal variables for [orm.fit()] and [lrm.fit()].  It was written because just using [factor()] creates slightly different distinct `y` levels on different hardware because [factor()] uses [unique()] which functions slightly differently on different systems when there are non-significant digits in floating point numbers.
#'
#' @title recode2integer
#' @param y a numeric, factor, or character vector with no `NA`s
#' @param precision number of places to the right of the decimal place to round `y` if `y` is numeric but not integer, for the purpose of finding the distinct values.  Real values rounding to the same values under `precision` are mapped to the same integer output `y`
#' @param ftable set to `FALSE` to suppress creation of `freq`
#'
#' @return a list with the following elements:
#'   * `y`: vector of integer-coded `y`
#'   * `ylevels`: vector of corresponding original `y` values, possibly rounded to `precision`.  This vector is numeric unless `y` is `factor` or character, in which case it is a character vector.
#'   * `freq`: frequency table of rounded or categorical `y`, with `names` attribute for the (possibly rounded) `y` levels of the frequencies
#'   * `median`: median `y` from original values if numeric, otherwise median of the new integer codes for `y`
#'   * `whichmedian`: the integer valued `y` that most closely corresponds to `median`; for an ordinal regression model this represents one plus the index of the intercept vector corresponding to `median`.
#' @export
#' @author Cole Beck
#' @md
#'
#' @examples
#' w <- function(y, precision=7) {
#'   v <- recode2integer(y, precision);
#'   print(v)
#'   print(table(y, ynew=v$y))
#' }
#' set.seed(1)
#' w(sample(1:3, 20, TRUE))
#' w(sample(letters[1:3], 20, TRUE))
#' y <- runif(20)
#' w(y)
#' w(y, precision=2)

recode2integer <- function(y, precision=7, ftable=TRUE) {

  # four scenarios for "y"
  # 1. y is numeric and contains decimals
  # 2. y is numeric and does not contain decimals
  # 3. y is factor/categorical
  # 4. y is something else (character)
  y_new    <- NULL
  ynumeric <- is.numeric(y)
  if(ynumeric) {
    # median of "y"
    mediany <- quantile(y, probs = 0.5, type = 7)
    # need precision if any fractional values
    needPrecision <- any(y %% 1 != 0)
    if(needPrecision) {
      ## scenario #1
      # when determining unique values of "y", round to avoid unpredictable behavior
      # this is better than `round(y, precision)`
      y_rnd       <- round(y       * 10^precision)
      mediany_rnd <- round(mediany * 10^precision)
      # distinct values of "y"
      yu <- sort(unique(y_rnd))
      # convert whole number back to decimal
      ylevels <- yu * 10^-precision
      # map "y" values from 1:n for `n` unique value
      y_new <- match(y_rnd, yu)
      # find the midpoint
      whichmedian <- which.min(abs(yu - mediany_rnd))
    } else {
      ## scenario #2
      yu      <- sort(unique(y))
      ylevels <- yu
      y_new   <- match(y, yu)
      whichmedian <- which.min(abs(yu - mediany))
    }
  }
  # For large n, as.factor is slow
  # if(!is.factor(y)) y <- as.factor(y)
  if(is.factor(y)) {
    ## scenario #3
    ylevels <- levels(y)
    y       <- as.integer(y)
  }
  else {
    if(length(y_new)) {
      # work already done if "y_new" is set
      y       <- y_new
    } else {
      ## scenario #4
      # if not done, map "y" values from 1:n for `n` unique value
      ylevels <- sort(unique(y))
      y       <- match(y, ylevels)
    }
  }
  if(! ynumeric) {
    yu      <- sort(unique(y))
    mediany <- quantile(y, probs = 0.5, type = 7)
    whichmedian <- which.min(abs(yu - mediany))
  }

  list(y=y, ylevels=ylevels,
       freq=if(ftable) structure(tabulate(y), names=ylevels),
       median=unname(mediany), whichmedian=whichmedian)
}

