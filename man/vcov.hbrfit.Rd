\name{vcov.hbrfit}
\alias{vcov.hbrfit}
\title{
Variance-covariance of regression estimates based on an HBR fit.
}
\description{ Variance-covariance of regression estimates based on an HBR fit.  }
\usage{
\method{vcov}{hbrfit}(object,...)
}
\arguments{
  \item{object}{object of class 'hbrfit'.  usually the result of a call to 'hbrfit'.}
  \item{\dots}{ additional arguments. currently unused.}
}
\value{
(p+1) x (p+1) variance covariance matrix where p is the ncol(x) and x is the 
design matrix.
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

\seealso{\code{\link{hbrfit}}}

\examples{
n<-40 ; p<-2
x<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)
vcov(hbrfit(y~x))
vcov(rfit(y~x))
}
