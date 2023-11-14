\name{psi}
\alias{psi}
\title{Psi function used in hbrfit}
\description{
Psi function used in the calculation of the weights for HBR estimates.
}
\usage{
psi(x)
}
\arguments{
  \item{x}{n x 1 vector}
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
Jeff Terpstra, Joe McKean, John Kloke 
}
\examples{
## The function is currently defined as
function (x) 
{
    x[x == -Inf] = -100
    x[x == Inf] = 100
    -1 * (x <= -1) + x * (-1 < x & x < 1) + 1 * (x >= 1)
  }
}
