\name{setPb}
\alias{setPb}
\title{Progress Bar for Simulations}
\description{
Depending on prevailing \code{options(showprogress=)} and availability
of the \code{tcltk} package, sets up a progress bar and creates a
function for simple updating of the bar as iterations progress.  Setting
\code{options(showprogressbar=FALSE)} or
\code{options(showprogressbar='none')} results in no progress being
shown.  Setting the option to something other than \code{"tk"} or
\code{"none"} results in the console being used to show the current
iteration number and intended number of iterations, the same as if
\code{tcltk} is not installed.  It is not recommended that the
\code{"tk"} be used for simulations requiring fewer than 10 seconds for
more than 100 iterations, as the time required to update the pop-up
window will be more than the time required to do the simulations.  This
problem can be solved by specifying, for example, \code{every=10} to
\code{setPb} or to the function created by \code{setPb}, or by using
\code{options(showevery=10)} before \code{setPb} is called.  If
\code{options(showprogress=)} is not specified, progress is shown in the
console with an iteration counter.
}
\usage{
setPb(n, type = c("Monte Carlo Simulation", "Bootstrap",
                  "Cross-Validation"),
         label, usetk = TRUE, onlytk=FALSE, every=1)
}
\arguments{
  \item{n}{maximum number of iterations}
  \item{type}{type of simulation.  Used for the progress bar title if
	\code{tcltk} is being used.}
  \item{label}{used to customize the bar label if present, overriding \code{type}}
  \item{usetk}{set to \code{FALSE} to override, acting as though the
	\code{tcltk} package were not installed}
  \item{onlytk}{set to \code{TRUE} to not write to the console even if
	\code{tcltk} is unavailable and \code{showprogressbar} is not
	\code{FALSE} or \code{"none"}}
  \item{every}{print a message for every \code{every} iterations}
}
\value{a function that should be called by the user once per iteration,
  specifying the iteration number as the sole argument}
\author{Frank Harrell}
\seealso{\code{\link[tcltk]{tkProgressBar}}, \code{\link[tcltk]{setTkProgressBar}}}
\examples{
\dontrun{
options(showprogress=TRUE)   # same as ='tk'
pb <- setPb(1000)
for(i in 1:1000) {
   pb(i)   # pb(i, every=10) to only show for multiples of 10
   # your calculations
  }
# Force rms functions to do simulations to not report progress
options(showprogress='none')
# For functions that do simulations to use the console instead of pop-up
# Even with tcltk is installed
options(showprogress='console')
pb <- setPb(1000, label='Random Sampling')
}
}
\keyword{utilities}
