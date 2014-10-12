\name{summary.hbrfit}
\alias{summary.hbrfit}
\title{ Summarize High Breakdown Rank-Based Linear Model Fits }
\description{
Provides a summary similar to the traditional least squares fit.
}
\usage{
\method{summary}{hbrfit}(object, ...)
}
\arguments{
  \item{object}{an object of class 'hbrfit', usually, a result of a call to 'hbrfit'}
  \item{\dots}{additonal arguments.  currently unused.}
}

\references{
Chang, W. McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999),
High breakdown rank-based regression,
\emph{Journal of the American Statistical Association}, 94, 205-219.

Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Terpstra, J. and McKean, J.W. (2005), Rank-based analyses of linear models using R,
\emph{Journal of Statistical Software}, 14(7).
}
\author{
Jeff Terpstra, Joe McKean,
John Kloke \email{kloke@biostat.wisc.edu}
}
\seealso{\code{\link{hbrfit}}}
\examples{
data(stars)
summary(hbrfit(light~temperature,data=stars))
}
