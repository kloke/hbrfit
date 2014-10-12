\name{stars}
\alias{stars}
\docType{data}
\title{Stars data.}
\description{
This dataset is drawn from an astronomy study on the star cluster CYG OB1 which contains 47 stars.
}
\usage{data(stars)}
\format{
  A data frame with 47 observations on the following 2 variables.
  \describe{
    \item{\code{temperature}}{The explanatory variable.  The logarithm of the temperature of the star.}
    \item{\code{light}}{The response variable. The logarithm of the light intensity of the star.}
  }
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Terpstra, J. and McKean, J.W. (2005), Rank-based analyses of linear models using R,
\emph{Journal of Statistical Software}, 14(7).
}
\examples{
data(stars)
plot(stars)
abline(hbrfit(light~temperature,data=stars))
}
\keyword{datasets}
