\name{hbrwts}
\alias{hbrwts}
\title{ Returns the weights for HBR regression.  }
\description{ Returns the weights for HBR regression.  }
\usage{
hbrwts(x, y, percent = 0.95, ehat0 = ltsreg(x, y)$residuals)
}
\arguments{
  \item{x}{n by p design matrix }
  \item{y}{n by 1 vector of responses }
  \item{percent}{ percentile of chi^2 dist to be used as a cut-off }
  \item{ehat0}{initial residuals }
}
\details{
Used internally in hbrfit.  For details see references below.
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
