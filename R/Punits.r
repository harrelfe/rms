#' Prepare units for Printing and Plotting
#'
#' Takes a character variable containing units of measurement for a variable.
#' If it has zero length, a `""` string is return.  Otherwise, any trailing `"s"` is
#' removed if the string is longer than one character, and depending on the arguments,
#' the string is changed to lower case, `"s"` is added, and the first character is
#' changed to upper case.
#'
#' @param u a single string containing units of measurement
#' @param lower if `TRUE` set string to all lower case
#' @param adds if `TRUE` add trailing `"s"`
#' @param upfirst if `TRUE` set first character to upper case
#' @param default default units if `u` is empty
#'
#' @returns a single character string
#' @seealso [Hmisc::units()]
#' @md
#'
#' @examples
#' \dontrun{
#' Punits('Years')
#' }
Punits <- function(u, lower=TRUE, adds=TRUE, upfirst=FALSE, default='') {
  if((! length(u) || u == '') && default != '') u <- default
  if(! length(u)) return('')
  
  if(lower) u <- tolower(u)
  if(nchar(u) > 1) {
    u <- sub('s$', '', u)
    if(adds) u <- paste0(u, 's')
    if(upfirst) u <- upFirst(u)
  }
  u
}
