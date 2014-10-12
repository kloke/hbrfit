\name{pairup}
\alias{pairup}
\title{Returns pairs of observations for use in hbrfit calculations.}
\description{
Internal function.  Returns pairs of observations.
Options for all pairwise, based on indices being less than to (i <= j), or indices 
being less than ( i < j).
}
\usage{
pairup(x, type = "less")
}
\arguments{
  \item{x}{n by p matrix }
  \item{type}{one of 'less', 'leq', 'neq'}
}
\details{
The calculation of all n^2 pairs is obtain and the requested subset is return.
The default being those with i < j.
The result is a matrix with 2p columns.
}
\author{
Jeff Terpstra, Joe McKean,
John Kloke \email{kloke@biostat.wisc.edu}
}
\examples{
## The function is currently defined as
function (x, type = "less") 
{
    x = as.matrix(x)
    n = dim(x)[1]
    i = rep(1:n, rep(n, n))
    j = rep(1:n, n)
    c1 = apply(x, 2, function(y) {
        rep(y, rep(length(y), length(y)))
    })
    c2 = apply(x, 2, function(y) {
        rep(y, length(y))
    })
    ans = cbind(c1, c2)
    switch(type, less = ans[(i < j), ], leq = ans[i <= j, ], 
        neq = ans)
  }
}
