\name{importedexported}
\alias{Surv}
\alias{ggplot}
\title{Exported Functions That Were Imported From Other Packages}
\description{
\code{Surv} and \code{ggplot} are imported from, respectively, the
\code{survival} and \code{ggplot2} packages and are exported from
\code{rms} so that the user does not have to attach these packages to do
simple things.
}

\usage{
Surv(time, time2, event,
     type = c("right", "left", "interval", "counting", "interval2", "mstate"),
     origin = 0)

ggplot(data = NULL, mapping = aes(), ..., environment =  parent.frame())
}

\arguments{
\item{time, time2, event, type, origin}{see \code{\link[survival]{Surv}}}
\item{data, mapping, ..., environment}{see \code{\link[ggplot2]{ggplot}}}
}
\value{see documentation in the original packages}
\seealso{
\code{\link[survival]{Surv}},
\code{\link[ggplot2]{ggplot}}
}
\examples{
\dontrun{
f <- psm(Surv(dtime, death) ~ x1 + x2 + sex + race, dist='gau')
ggplot(Predict(f))
}
}
