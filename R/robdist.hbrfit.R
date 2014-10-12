# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
robdist.hbrfit <-
function (x) 
{
    x.o <- x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    qn <- floor((n + p + 1)/2)
    if (qn < p + 1) 
        stop(paste("quantile must be at least", p + 1))
    divisor <- apply(x, 2, IQR)
    if (any(divisor == 0)) 
        stop("mycov.rob:  at least one column has IQR 0")
    x <- x/rep(divisor, rep(n, p))
    if (p > 1) {
        best <- covMcd(x)$best
        if (!length(best)) {
            stop("x is probably collinear")
        }
        else {
            means <- colMeans(x[best, , drop = FALSE])
            rcov <- var(x[best, , drop = FALSE]) * (1 + 15/(n - 
                p))^2
        }
    }
    else {
        z <- covMcd(x)
        means <- z$center
        rcov <- z$cov
    }
    dist <- mahalanobis(x, means, rcov)
    cut <- qchisq(0.975, p) * quantile(dist, qn/n)/qchisq(qn/n, 
        p)
    center = colMeans(x[dist < cut, , drop = FALSE]) * divisor
    cov <- divisor * var(x[dist < cut, , drop = FALSE]) * rep(divisor, 
        rep(p, p))
    mahalanobis(x.o, center, cov)
}
