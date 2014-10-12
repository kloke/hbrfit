\name{hbrfit}
\alias{hbrfit}
\title{
High Breakdown Rank (HBR) Estimates
}
\description{
High breakdown rank (HBR) estimates are robust to outliers in 
both X & Y spaces.  They are based on a weighted Wilcoxon pseudo-norm.
Data points which are outliers in both X & Y space are downweighted.
HBR estimates achieve 50% breakdown point.
}
\usage{
hbrfit(formula, data, subset, symmetric = FALSE,...)
}
\arguments{
  \item{formula}{an object of class formula}
  \item{data}{an optional data frame}
  \item{subset}{an optional argument specifying the subset of observations to be used}
  \item{symmetric}{logical.  If 'FALSE' uses median of residuals as estimate of intercept}
  \item{\dots}{additional arguments.  Currently unused.}
}
\details{
The HBR pseudo-norm is ||u|| = sum_{i < j} b_{ij} |u_i - u_j|.
The weights (b_{ij}) are chosen based on robust measures of distance.
Data points with large residuals (based on a initial LTS fit) and outling design points
(based on a robust measure of Mahalanobis distance) are down weighted (b_{ij} < 1).
If all b_{ij} = 1 then the HBR pseudo-norm is the Wilcoxon.
HBR estimates for linear models were developed by Chang, et. al. (1999).  See also Section 3.12 of Hettmansperger and McKean (2011).
}
\value{
\item{coefficients}{estimated regression coefficents with intercept}
\item{residuals}{the residuals, i.e. y-yhat}
\item{fitted.values}{ yhat = x betahat}
\item{weights}{estimated weights. the b_{ij}}
\item{x}{original design matrix}
\item{y}{original response vector}
\item{tauhat}{estimated value of the scale parameter tau}
\item{taushat}{estimated value of the scale parameter tau_s}
\item{betahat}{estimated regression coefficents}
\item{qrx1}{qrd of the design matrix with a column of ones prepended}
\item{call}{Call to the function}
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

\seealso{\code{\link{summary.hbrfit}}}

\examples{
data(stars)
hbrfit(light~temperature,data=stars)
}
\keyword{nonparametric}
\keyword{robust}
\keyword{regression}
