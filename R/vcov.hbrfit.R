# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
vcov.hbrfit <-
function (object,...) 
{
	fit<-object
    x <- as.matrix(with(fit, qr.Q(qrx1)[, 2:qrx1$rank]))
    bmat = as.matrix(fit$weights)
    res = as.vector(fit$residuals)
    n = dim(x)[1]
    p = dim(x)[2]
    tau <- fit$tauhat
    diag(bmat) = rep(0, n)
    w = -1 * bmat
    diag(w) = bmat %*% as.matrix(rep(1, n))
    w = (1/(sqrt(12) * tau)) * w
    cmat = (1/n^2) * t(x) %*% w %*% x
    cinv = solve(cmat)
    u = (1/n) * (bmat - diag(c(bmat %*% cbind(rep(1, n))))) %*% 
        x
    u = u * (1 - 2 * rank(res)/n)
    vmat = var(u)
    varcovbetahatstar = (1/(4 * n)) * cinv %*% vmat %*% cinv
    Q <- qr.Q(fit$qrx1)
    q1 <- Q[, 1]
    Q2 <- Q[, 2:fit$qrx1$rank]
    xxpxi <- fit$x %*% chol2inv(chol(crossprod(fit$x)))
    A1 <- crossprod(q1, xxpxi)
    A2 <- crossprod(Q2, xxpxi)
    varcov <- fit$taushat^2 * crossprod(A1) + t(A2) %*% varcovbetahatstar %*% 
        A2
    attr(varcov, "names") = NULL
    varcov
}
