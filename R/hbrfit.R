# Based on weighted Wilcoxon code developed by Terpstra & McKean (2005) 
# (http://www.jstatsoft.org/v14/i07) under GPL 2.
hbrfit <-
function (formula, data, subset, symmetric = FALSE,...) 
{
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    x.o <- x <- model.matrix(attr(mf, "terms"), data = mf)
    y.o <- y <- model.response(mf)
    x <- as.matrix(x[, colnames(x) != "(Intercept)"])
    x1 <- cbind(rep(1, nrow(x)), x)
    qrx1 <- qr(x1)
    Q <- as.matrix(qr.Q(qrx1))
    q1 <- Q[, 1]
    x <- as.matrix(Q[, 2:qrx1$rank])
    n <- nrow(x)
    p <- ncol(x)
    bij <- hbrwts(x, y)
    ypairs = pairup(y)
    yi = ypairs[, 1]
    yj = ypairs[, 2]
    xpairs = pairup(x)
    xi = xpairs[, 1:p]
    xj = xpairs[, (p + 1):(2 * p)]
    ystar = bij * (yi - yj)
    xstar = bij * (xi - xj)
    est <- rq(ystar ~ xstar - 1, method = ifelse(length(ystar) > 
        5000 | p > 20, "fnb", "br"))$coefficients
    int = median(y - (x %*% as.matrix(est)))
    resid = as.vector(y - int - (x %*% as.matrix(est)))
    yhat <- y - resid
    bhat <- lsfit(x.o, yhat, intercept = F)$coefficients
    wts = matrix(0, nrow = n, ncol = n)
    index = pairup(1:n)
    wts[index] = bij
    wts[index[, 2:1]] = bij
    tauhat <- gettauF0(resid, p)
    if (symmetric) {
        taushat <- tauhat
    }
    else {
        taushat <- taustar(resid, p)
    }
    ans = list(coefficients = bhat, residuals = resid, fitted.values = y - 
        resid, weights = wts, y = y.o, x = x.o, tauhat = tauhat, 
        taushat = taushat, betahat = bhat, qrx1 = qrx1)
    ans$call <- call
    class(ans) <- c("hbrfit", "rfit")
    ans
}
