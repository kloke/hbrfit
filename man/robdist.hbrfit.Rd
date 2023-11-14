\name{robdist.hbrfit}
\alias{robdist.hbrfit}
\title{Robust distances.  Internal function for hbrfit}
\description{
Returns the robust measures of distance used in the calculation of the weights.
Based on a robust measure of Mahalanobis distance.  Center and variance 
are estimated based on a high breakdown point Minimum Covariance Determinate 
estimate.
}
\usage{ robdist.hbrfit(x) }
\arguments{ \item{x}{n by p matrix.  usually the design matrix.} }
\references{
Chang, W. McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999),
High breakdown rank-based regression,
\emph{Journal of the American Statistical Association}, 94, 205-219.

Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Terpstra, J. and McKean, J.W. (2005), Rank-based analyses of linear models using R,
\emph{Journal of Statistical Software}, 14(7).
}
\author{
Jeff Terpstra, Joe McKean, John Kloke 
}

