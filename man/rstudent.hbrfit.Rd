\name{rstudent.hbrfit}
\alias{rstudent.hbrfit}
\title{Studentized Residuals for HBR Regression}
\description{
Returns the Studentized residuals based on high breakdown rank-based estimation.
}
\usage{
\method{rstudent}{hbrfit}(model, ...)
}
\arguments{
  \item{model}{an object of class hbrfit}
  \item{\dots}{additonal arguments.  currently not used.}
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
x<-runif(47)
y<-rcauchy(47)
qqnorm(rstudent(fit<-rfit(y~x)))
plot(x,rstudent(fit)) ; abline(h=c(-2,2))
}
