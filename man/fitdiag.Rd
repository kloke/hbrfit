\name{fitdiag}
\alias{fitdiag}
\title{ Robust Fit Diagnostics }
\description{
Diagnostics to compare two fits of a linear model.   The comparison (tdbeta) is a measure of total difference between the fits while the comparison (cfit) is measure of the difference between the fitted (predicted) values.   Comparisons between any two of the following fits can be made: Wilcoxon (WIL), least squares (LS), HBR (HBR), and least trim squares.
}
\usage{
fitdiag(x, y, est = c("WIL", "HBR"), ...)
}
\arguments{
  \item{x}{design matrix.}
  \item{y}{response vector.}
  \item{est}{character vector of length 2.  estimates to use. default is est = c("WIL", "HBR"))}
  \item{\dots}{additional arguments passed to rfit.}
}
\details{
The diagnostic tdbeta uses the benchmark bmtd.  Values exceeding bmtd indicate the fits differ.  The diagnostic cfit uses the benchmark bmcf; for each case_i, cfit_i exceeding bmcf in absolute value indicates the ith fitted values differ.  These are only benchmarks, not critical values.   A useful plot is cfit versus case.  The user should look for large gaps in this plot.  
}
\value{
	\item{tdbeta}{Total difference.}
	\item{cfit}{Casewise differences in fits.}
	\item{bmtd}{Benchmark for tdbeta.}
	\item{bmcf}{Benchmark for cfit.}
	\item{est}{The two fits being compared.}
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{
Jeff Terpstra, Joe McKean, John Kloke 
}
\examples{
# Case of no outliers
set.seed(35443092)
x<-rnorm(20);y<-rnorm(20)
dwl <- fitdiag(x,y,est=c("WIL","LS"))
tdbetas <- round(c(dwl$tdbeta,dwl$bmtd),digits=2)
tdbetas
Case<-1:length(y)
plot(dwl$cfit~Case,ylab='Difference in fits',ylim=c(-dwl$bmcf,dwl$bmcf))
title('Simulated Normal Data')
abline(h=c(-dwl$bmcf,dwl$bmcf))  # benchmark bars
# Case where there is a gap
dwilhbr = fitdiag(stars$temperature,stars$light,est=c("WIL","HBR"))
tdbetas <- round(c(dwilhbr$tdbeta,dwilhbr$bmtd),digits=2)
tdbetas
Case <- 1:length(stars$temperature)
plot(dwilhbr$cfit~Case,ylab='Difference in fits')
title('Stars Data with Six Outliers')
}
\keyword{robust}

