\name{diagplot}
\alias{diagplot}
\title{Diagnostic plots based on a linear model fit.}
\description{
Diagnostic plots based on a linear model fit.
Used for it's side effect of creating a 2x2 grid with the 
following plots: raw residuals vs fits; histogram of raw residuals;
case plot of the Studentized residuals; normal probability plot of the raw residuals.
}
\usage{diagplot(fit, ...)}
\arguments{
  \item{fit}{a linear model fit.}
  \item{\dots}{additional arguments.  currently unused.}
}
\details{
Designed for use with hbrfit and rfit, but should work with any linear model fit such that
the functions fitted.values(fit), residuals(fit), rstudent(fit) return the appropriate 
vectors.
}
\references{
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
# R estimates are robust to outliers in Y-space
x<-runif(47)
y<-rcauchy(47)
fit<-rfit(y~x)
diagplot(fit)

# HBR estimates are robust to outliers in both X- and Y-space
x<-rcauchy(47)
y<-rcauchy(47)
diagplot(hbrfit(y~x))
}
