\name{hbrfit-package}
\alias{hbrfit-package}
\docType{package}
\title{High Breakdown Rank-based (HBR) fitting and inference for linear regression models}
\description{
High breakdown rank-based (HBR) estimates are robust to outliers in
both X & Y spaces.  They are based on a weighted Wilcoxon pseudo-norm.
Data points which are outliers in both X & Y space are downweighted.
HBR estimates achieve 50% breakdown point.

This package is based on the weighted Wilcoxon (ww) code developed by
Terpstra and McKean (2005) under GPL.
}
\details{
\tabular{ll}{
Package: \tab hbrfit\cr
Type: \tab Package\cr
Version: \tab 0.2.6\cr
Date: \tab 2023-11-14\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Jeff Terpstra, Joe McKean, John Kloke 
}
\references{
Chang, W. McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999),
High breakdown rank-based regression,
\emph{Journal of the American Statistical Association}, 94, 205-219.

Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Terpstra, J. and McKean, J.W. (2005), Rank-based analyses of linear models using R,
\emph{Journal of Statistical Software}, 14(7).
}
\keyword{nonparametric}
\keyword{robust}
\keyword{regression}
\examples{
data(stars)
plot(stars)
fit<-hbrfit(light~temperature,data=stars)
abline(fit)
summary(fit)
}
