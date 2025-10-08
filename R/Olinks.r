#' Likehood-Based Statistics for Other Links for orm Fits
#'
#' @param object an object created by `orm` with `x=TRUE, y=TRUE`
#' @param links a vector of links to consider other than the one used to get `object`
#' @param dec number of digits to the right of the decimal place to round statistics to
#' @param gradtol tolerance for convergence on the absolute gradient; see `lrm.fit` and `orm.fit`.
#'
#' @returns data frame.  The `R2` column is from the last adjusted \eqn{R^2} computed by `orm`,
#' which adjustes for the effective sample size and the number of betas.
#' @export
#' @md
#' @author Frank Harrell
#'
#' @examples
#' \dontrun{
#' f <- orm(y ~ x1 + x2, family='loglog', x=TRUE, y=TRUE)
#' Olinks(f)
#' }
Olinks <- function(object, links=c('logistic', 'probit', 'loglog', 'cloglog'), dec=3, gradtol=0.001) {
  if(! inherits(object, 'orm')) stop('object must an orm object')
  if(! all(c('x', 'y') %in% names(object))) stop('must run orm with x=TRUE, y=TRUE')
  fam    <- object$family
  links  <- unique(c(fam, links))
  p      <- length(coef(object))
  fitter <- quickRefit(object, storevals=TRUE, compstats=TRUE, gradtol=gradtol, what='fitter')
  R      <- NULL
  for(fm in links) {
    f    <- if(fm ==fam) object else fitter(family=fm)
    dev  <- deviance(f)
    st   <- f$stats
    r2   <- st[grep('^R2', names(st))]
    last <- length(r2)
    r2   <- unname(r2[last])
    aic  <- dev[2] + 2 * p
    d    <- data.frame(link=fm, null.deviance=dev[1], deviance=dev[2],
                       AIC=aic, LR=st['Model L.R.'], R2=r2)
    R   <- rbind(R, d)
  }
  row.names(R) <- NULL
  for(i in 2 : 6) R[[i]] <- round(R[[i]], dec)
  R
}
